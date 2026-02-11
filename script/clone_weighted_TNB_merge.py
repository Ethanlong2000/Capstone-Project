#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build long-format TNB table for B2:
- Supports both IC50-based and %Rank-based neoantigen definitions
- Counting units: mutation / peptide / hla_peptide
- Outputs: count TNB + weighted TNB, each with total and clonal (VAF cutoffs)

Output long-format columns:
sample, unit, metric_family, strategy, binder_metric, clonal_cutoff, value, n_total_unique
"""

import argparse
import glob
import os
import re
from typing import Dict, List, Optional, Tuple
import math

import pandas as pd


# -----------------------------
# Utilities
# -----------------------------
def _safe_float(x):
    try:
        if x is None:
            return float("nan")
        x = str(x).strip()
        if x == "" or x.lower() in {"na", "nan", "none"}:
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def normalize_hla(allele: str) -> str:
    if allele is None:
        return ""
    a = str(allele).strip()
    if a == "" or a.lower() in {"na", "nan", "none"}:
        return ""
    a = a.upper().replace(" ", "")
    a = a[4:] if a.startswith("HLA-") else a

    m = re.match(r"^([ABC])\*?(\d{2,3})(:?)(\d{2,3})$", a)
    if m:
        locus = m.group(1)
        f1 = m.group(2)
        f2 = m.group(4)
        return f"HLA-{locus}*{f1}:{f2}"

    m2 = re.match(r"^([ABC])\*(\d{2,3}):(\d{2,3})$", a)
    if m2:
        return f"HLA-{m2.group(1)}*{m2.group(2)}:{m2.group(3)}"

    if re.match(r"^[ABC].*", a):
        return "HLA-" + a
    return a


def get_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def build_mutation_id(df: pd.DataFrame) -> pd.Series:
    if "Mutation" in df.columns:
        return df["Mutation"].astype(str).fillna("").str.strip()

    needed = ["Chromosome", "Start", "Reference", "Variant"]
    for c in needed:
        if c not in df.columns:
            return df.index.astype(str)

    return (
        df["Chromosome"].astype(str).fillna("").str.strip()
        + ":"
        + df["Start"].astype(str).fillna("").str.strip()
        + ":"
        + df["Reference"].astype(str).fillna("").str.strip()
        + ">"
        + df["Variant"].astype(str).fillna("").str.strip()
    )


def collect_files(input_path: str, pattern: str = "*.tsv") -> List[str]:
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        return sorted(glob.glob(os.path.join(input_path, pattern)))
    raise FileNotFoundError(f"未找到路径: {input_path}")


# -----------------------------
# Dedup keys and "best binder" collapsing
# -----------------------------
def unit_keys(unit: str) -> List[str]:
    if unit == "mutation":
        return ["mutation_id"]
    if unit == "peptide":
        return ["mutation_id", "mt_pep"]
    if unit == "hla_peptide":
        return ["mt_pep", "hla_norm"]
    raise ValueError(f"Unknown unit: {unit}")


def collapse_best_per_unit(
    df: pd.DataFrame,
    unit: str,
    prefer: str = "rank",
) -> pd.DataFrame:
    """
    Important for mutation/peptide units:
    - same mutation generates many peptides & HLA rows
    - we collapse to ONE representative per unit (best binder) before counting/weighting
    For hla_peptide unit, rows are already at the unit-level -> no collapse.

    prefer:
      - "rank": choose row with minimal mt_rank, tie-breaker minimal mt_ic50
      - "ic50": choose row with minimal mt_ic50, tie-breaker minimal mt_rank
    """
    keys = unit_keys(unit)
    if unit == "hla_peptide":
        return df.drop_duplicates(subset=keys, keep="first")

    work = df.copy()
    # sort so the "best" appears first within each group
    if prefer == "rank":
        work = work.sort_values(
            by=["mt_rank", "mt_ic50"],
            ascending=[True, True],
            na_position="last",
            kind="mergesort",
        )
    else:
        work = work.sort_values(
            by=["mt_ic50", "mt_rank"],
            ascending=[True, True],
            na_position="last",
            kind="mergesort",
        )
    return work.drop_duplicates(subset=keys, keep="first")


# -----------------------------
# Weighted scoring (stable, not exploding)
# -----------------------------
def weight_bind_from_rank(rank: float) -> float:
    # w_bind = max(0, -log10(rank)) ; rank=0.01 -> 2 ; rank=1 -> 0
    if pd.isna(rank) or rank <= 0:
        return 0.0
    import math
    return max(0.0, -math.log10(rank))



def weight_bind_from_ic50(ic50: float, ref: float = 500.0) -> float:
    # monotonic: stronger binder -> larger weight
    # use -log10(IC50) around nM scale, clamp negatives
    if pd.isna(ic50) or ic50 <= 0:
        return 0.0
    return max(0.0, math.log10(ref / ic50))

def weight_expr_from_tpm(tpm: float) -> float:
    # log2(TPM+1)
    if pd.isna(tpm) or tpm < 0:
        return 0.0
    import math
    return math.log2(tpm + 1.0)


def weight_clonal_from_vaf(vaf: float) -> float:
    # min(1, VAF/0.5)
    if pd.isna(vaf) or vaf < 0:
        return 0.0
    return min(1.0, vaf / 0.5)


def compute_weighted_score_row(
    binder_metric: str,
    mt_rank: float,
    mt_ic50: float,
    tpm: float,
    vaf: float,
) -> float:
    if binder_metric == "rank":
        w_bind = weight_bind_from_rank(mt_rank)
    else:
        w_bind = weight_bind_from_ic50(mt_ic50)
    w_expr = weight_expr_from_tpm(tpm)
    w_clonal = weight_clonal_from_vaf(vaf)
    return w_bind * w_expr * w_clonal


# -----------------------------
# Strategy definitions
# -----------------------------
def build_strategies() -> Dict[str, Dict]:
    """
    Keep naming close to your existing S1-S4:
    IC50:
      S1: mt_ic50 < 500
      S2: mt_ic50 < 500 & TPM > 1
      S3: mt_ic50 < 50  & TPM > 5
      S4: mt_ic50 < 50  & TPM > 5 & VAF > 0.1 & wt_ic50 > 1000
    Rank:
      R1: mt_rank < 2
      R2: mt_rank < 2 & TPM > 1
      R3: mt_rank < 0.5 & TPM > 5
      R4: mt_rank < 0.5 & TPM > 5 & VAF > 0.1 & wt_rank > 20   (fallback if wt_rank exists)
          if wt_rank missing, will fallback to wt_ic50 > 1000 if wt_ic50 exists
    """
    return {
        "S1": {"binder_metric": "ic50", "mt_ic50_lt": 500, "tpm_gt": None, "vaf_gt": None, "wt_filter": None},
        "S2": {"binder_metric": "ic50", "mt_ic50_lt": 500, "tpm_gt": 1,    "vaf_gt": None, "wt_filter": None},
        "S3": {"binder_metric": "ic50", "mt_ic50_lt": 50,  "tpm_gt": 5,    "vaf_gt": None, "wt_filter": None},
        "S4": {"binder_metric": "ic50", "mt_ic50_lt": 50,  "tpm_gt": 5,    "vaf_gt": 0.1,  "wt_filter": ("wt_ic50_gt", 1000)},
        "R1": {"binder_metric": "rank", "mt_rank_lt": 2,   "tpm_gt": None, "vaf_gt": None, "wt_filter": None},
        "R2": {"binder_metric": "rank", "mt_rank_lt": 2,   "tpm_gt": 1,    "vaf_gt": None, "wt_filter": None},
        "R3": {"binder_metric": "rank", "mt_rank_lt": 0.5, "tpm_gt": 5,    "vaf_gt": None, "wt_filter": None},
        "R4": {"binder_metric": "rank", "mt_rank_lt": 0.5, "tpm_gt": 5,    "vaf_gt": 0.1,  "wt_filter": ("wt_rank_gt", 10)},
    }   


def apply_strategy_mask(df: pd.DataFrame, strat: Dict) -> pd.Series:
    m = pd.Series(True, index=df.index)

    if strat["binder_metric"] == "ic50":
        thr = strat.get("mt_ic50_lt", None)
        if thr is not None:
            m &= df["mt_ic50"] < thr
    else:
        thr = strat.get("mt_rank_lt", None)
        if thr is not None:
            m &= df["mt_rank"] < thr

    tpm_gt = strat.get("tpm_gt", None)
    if tpm_gt is not None:
        m &= df["tpm"] > tpm_gt

    vaf_gt = strat.get("vaf_gt", None)
    if vaf_gt is not None:
        m &= df["vaf"] > vaf_gt

    wt_filter = strat.get("wt_filter", None)
    if wt_filter is not None:
        k, v = wt_filter
        if k == "wt_rank_gt":
            # if wt_rank missing, fallback to wt_ic50 if present
            if df["wt_rank"].notna().any():
                m &= df["wt_rank"] > float(v)
            else:
                m &= df["wt_ic50"] > 1000
        elif k == "wt_ic50_gt":
            m &= df["wt_ic50"] > float(v)

    return m


# -----------------------------
# Core: process one file -> long rows
# -----------------------------
def process_file_to_long(
    path: str,
    units: List[str],
    vaf_cutoffs: List[float],
    sample_fullname: bool = False,
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    base = os.path.basename(path)
    name_no_ext = base.rsplit(".", 1)[0]
    sample = name_no_ext if sample_fullname else name_no_ext.split("_")[0]

    # Identify key columns
    mt_pep_col = get_col(df, ["MT Epitope Seq", "MT Epitope", "MT_peptide", "MT_PEPTIDE"])
    hla_col = get_col(df, ["HLA Allele", "HLA", "HLA_Allele", "Allele"])

    mt_ic50_col = get_col(
        df,
        [
            "NetMHCpan MT IC50 Score",
            "MT IC50",
            "mt_ic50",
            "MT_IC50",
        ],
    )
    wt_ic50_col = get_col(
        df,
        [
            "NetMHCpan WT IC50 Score",
            "WT IC50",
            "wt_ic50",
            "WT_IC50",
        ],
    )

    # IMPORTANT: %Rank columns vary a lot; add candidates here as needed
    mt_rank_col = get_col(
        df,
        [
            "NetMHCpan MT Percentile",
            "NetMHCpan MT %Rank",
            "NetMHCpan MT Rank",
            "MT Percentile",
            "MT %Rank",
            "MT Rank",
            "mt_rank",
            "MT_RANK",
        ],
    )
    wt_rank_col = get_col(
        df,
        [
            "NetMHCpan WT Percentile",
            "NetMHCpan WT %Rank",
            "NetMHCpan WT Rank",
            "WT Percentile",
            "WT %Rank",
            "WT Rank",
            "wt_rank",
            "WT_RANK",
        ],
    )

    vaf_col = get_col(df, ["Tumor DNA VAF", "VAF", "tumor_vaf", "TumorVAF"])
    tpm_col = get_col(df, ["TPM", "tpm", "Gene TPM", "RNA TPM"])

    if mt_pep_col is None:
        raise ValueError("Missing MT peptide column (e.g., 'MT Epitope Seq').")
    if vaf_col is None:
        raise ValueError("Missing VAF column (e.g., 'Tumor DNA VAF').")
    if mt_ic50_col is None and mt_rank_col is None:
        raise ValueError("Need at least one binder metric column: MT IC50 or MT %Rank.")

    work = df.copy()
    work["sample"] = sample
    work["file"] = base
    work["mutation_id"] = build_mutation_id(work)

    work["mt_pep"] = work[mt_pep_col].astype(str).fillna("").str.strip()
    work["hla_norm"] = work[hla_col].apply(normalize_hla) if hla_col else ""

    work["mt_ic50"] = work[mt_ic50_col].apply(_safe_float) if mt_ic50_col else float("nan")
    work["wt_ic50"] = work[wt_ic50_col].apply(_safe_float) if wt_ic50_col else float("nan")

    work["mt_rank"] = work[mt_rank_col].apply(_safe_float) if mt_rank_col else float("nan")
    work["wt_rank"] = work[wt_rank_col].apply(_safe_float) if wt_rank_col else float("nan")

    work["vaf"] = work[vaf_col].apply(_safe_float)
    work["tpm"] = work[tpm_col].apply(_safe_float) if tpm_col else float("nan")

    strategies = build_strategies()

    long_rows = []

    for unit in units:
        # total_unique (after collapsing best per unit) can be useful denominator
        # we compute it twice: prefer rank and prefer ic50 (minor difference). Use prefer=rank for reporting.
        total_unique_rank = len(collapse_best_per_unit(work, unit, prefer="rank"))
        total_unique_ic50 = len(collapse_best_per_unit(work, unit, prefer="ic50"))

        for strat_name, strat in strategies.items():
            binder_metric = strat["binder_metric"]

            # if strategy needs rank but mt_rank missing -> skip
            if binder_metric == "rank" and not work["mt_rank"].notna().any():
                continue
            if binder_metric == "ic50" and not work["mt_ic50"].notna().any():
                continue

            # collapse to best per unit BEFORE applying strategy, consistent and avoids peptide inflation
            collapsed = collapse_best_per_unit(work, unit, prefer=binder_metric)

            mask = apply_strategy_mask(collapsed, strat)
            sel = collapsed[mask].copy()

            # count: total (no clonal cutoff)
            n_sel = len(sel)
            long_rows.append({
                "sample": sample,
                "unit": unit,
                "metric_family": "count",
                "binder_metric": binder_metric,
                "strategy": strat_name,
                "clonal_cutoff": "all",
                "value": float(n_sel),
                "n_total_unique": float(total_unique_rank if binder_metric == "rank" else total_unique_ic50),
            })

            # weighted: total (requires TPM; if TPM missing -> still compute with expr=0 => all 0, so skip)
            if sel["tpm"].notna().any():
                sel["w"] = sel.apply(
                    lambda r: compute_weighted_score_row(
                        binder_metric=binder_metric,
                        mt_rank=r["mt_rank"],
                        mt_ic50=r["mt_ic50"],
                        tpm=r["tpm"],
                        vaf=r["vaf"],
                    ),
                    axis=1
                )
                long_rows.append({
                    "sample": sample,
                    "unit": unit,
                    "metric_family": "weighted",
                    "binder_metric": binder_metric,
                    "strategy": strat_name,
                    "clonal_cutoff": "all",
                    "value": float(sel["w"].sum()),
                    "n_total_unique": float(total_unique_rank if binder_metric == "rank" else total_unique_ic50),
                })

            # clonal versions across VAF cutoffs
            for t in vaf_cutoffs:
                sel_c = sel[sel["vaf"] >= float(t)]
                long_rows.append({
                    "sample": sample,
                    "unit": unit,
                    "metric_family": "count",
                    "binder_metric": binder_metric,
                    "strategy": strat_name,
                    "clonal_cutoff": f"vaf>={t}",
                    "value": float(len(sel_c)),
                    "n_total_unique": float(total_unique_rank if binder_metric == "rank" else total_unique_ic50),
                })

                if sel_c.shape[0] > 0 and sel_c["tpm"].notna().any():
                    # reuse weight
                    if "w" not in sel_c.columns:
                        sel_c = sel_c.copy()
                        sel_c["w"] = sel_c.apply(
                            lambda r: compute_weighted_score_row(
                                binder_metric=binder_metric,
                                mt_rank=r["mt_rank"],
                                mt_ic50=r["mt_ic50"],
                                tpm=r["tpm"],
                                vaf=r["vaf"],
                            ),
                            axis=1
                        )
                    long_rows.append({
                        "sample": sample,
                        "unit": unit,
                        "metric_family": "weighted",
                        "binder_metric": binder_metric,
                        "strategy": strat_name,
                        "clonal_cutoff": f"vaf>={t}",
                        "value": float(sel_c["w"].sum()),
                        "n_total_unique": float(total_unique_rank if binder_metric == "rank" else total_unique_ic50),
                    })

    return pd.DataFrame(long_rows)


def main():
    parser = argparse.ArgumentParser(description="Generate long-format TNB table (count + weighted + clonal) for B2.")
    parser.add_argument("input", help="QC 后的 TSV 文件或包含 TSV 的目录")
    parser.add_argument("--pattern", default="*.tsv", help="输入为目录时匹配模式，默认 *.tsv")
    parser.add_argument(
        "--units",
        default="mutation,peptide,hla_peptide",
        help="计数单位，逗号分隔：mutation,peptide,hla_peptide（默认全算）",
    )
    parser.add_argument(
        "--vaf-cutoffs",
        default="0.2,0.25,0.3,0.35",
        help="clonal VAF阈值列表，逗号分隔（默认 0.2,0.25,0.3,0.35）",
    )
    parser.add_argument(
        "--sample-fullname",
        action="store_true",
        help="使用去扩展名的完整文件名作为sample（默认取下划线前第一段）。",
    )
    parser.add_argument(
        "--output",
        help="输出CSV路径。默认写到输入目录/文件所在目录：TNB_long_B2.csv",
    )
    args = parser.parse_args()

    files = collect_files(args.input, args.pattern)
    if not files:
        print("未找到匹配的文件")
        return

    units = [u.strip() for u in args.units.split(",") if u.strip()]
    for u in units:
        if u not in {"mutation", "peptide", "hla_peptide"}:
            raise ValueError(f"--units 包含未知单位: {u}")

    vaf_cutoffs = [float(x.strip()) for x in args.vaf_cutoffs.split(",") if x.strip()]

    all_long = []
    for f in files:
        try:
            out = process_file_to_long(
                path=f,
                units=units,
                vaf_cutoffs=vaf_cutoffs,
                sample_fullname=args.sample_fullname,
            )
            all_long.append(out)
        except Exception as e:
            print(f"处理文件失败 {f}: {e}")

    if not all_long:
        print("没有成功生成任何结果")
        return

    long_df = pd.concat(all_long, ignore_index=True)

    # output path
    if args.output:
        out_path = args.output
    else:
        base_dir = os.path.dirname(args.input) if os.path.isfile(args.input) else args.input
        out_path = os.path.join(base_dir, "TNB_long_B2.csv")

    long_df.to_csv(out_path, index=False)
    print(f"已保存 long-format: {out_path}")
    print(long_df.head(20))


if __name__ == "__main__":
    main()

# python3 /work/longyh/BY/script/clone_weighted_TNB_merge.py qc_dir/ --pattern "*.tsv" \
#   --units mutation,peptide,hla_peptide \
#   --vaf-cutoffs 0.2,0.25,0.3,0.35
