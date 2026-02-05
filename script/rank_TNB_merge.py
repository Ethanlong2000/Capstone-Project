#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute TNB under multiple strategies using NetMHCpan %Rank (Percentile),
and output in long-format for easier plotting.

Counting units (dedup keys):
- mutation: unique mutation_id
- peptide: unique mutation_id + MT peptide
- hla_peptide: unique MT peptide + HLA allele

Strategies (Rank-based):
S1_rank: mt_rank < 2.0
S2_rank: mt_rank < 2.0 and TPM > 1
S3_rank: mt_rank < 0.5 and TPM > 5
S4_rank: mt_rank < 0.5 and TPM > 5 and VAF > 0.1 and wt_rank > 2.0

Notes:
- If TPM column is missing, only S1_rank will be computed.
- Default units: mutation (to keep the main analysis clean); you can pass --units to add more.
"""

import argparse
import glob
import os
import re
from typing import Dict, List, Optional

import pandas as pd


# -----------------------------
# Utilities
# -----------------------------
def _safe_float(x):
    try:
        if x is None:
            return float("nan")
        s = str(x).strip()
        if s == "" or s.lower() in {"na", "nan", "none"}:
            return float("nan")
        return float(s)
    except Exception:
        return float("nan")


def normalize_hla(allele: str) -> str:
    """
    Normalize common HLA-I allele formats to: HLA-A*02:01-like.
    Handles: HLA-A02:01 / HLA-A*02:01 / A*02:01 / A02:01 / a0201 / etc.
    If cannot parse, returns stripped original.
    """
    if allele is None:
        return ""
    a = str(allele).strip()
    if a == "" or a.lower() in {"na", "nan", "none"}:
        return ""

    a = a.upper().replace(" ", "")
    a = a[4:] if a.startswith("HLA-") else a

    # compact forms like A0201, B0702, C0602 or A02:01
    m = re.match(r"^([ABC])\*?(\d{2,3})(:?)(\d{2,3})$", a)
    if m:
        locus = m.group(1)
        f1 = m.group(2)
        f2 = m.group(4)
        return f"HLA-{locus}*{f1}:{f2}"

    # already like A*02:01
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
    """
    Prefer existing 'Mutation' column; else use Chromosome+Start+Reference+Variant.
    """
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


# -----------------------------
# Dedup keys by counting unit
# -----------------------------
def dedup_by_unit(df: pd.DataFrame, unit: str) -> pd.DataFrame:
    if unit == "mutation":
        return df.drop_duplicates(subset=["mutation_id"], keep="first")
    if unit == "peptide":
        return df.drop_duplicates(subset=["mutation_id", "mt_pep"], keep="first")
    if unit == "hla_peptide":
        return df.drop_duplicates(subset=["mt_pep", "hla_norm"], keep="first")
    raise ValueError(f"Unknown unit: {unit}")


# -----------------------------
# Core computation (Rank-based)
# -----------------------------
def compute_tnb_rank_counts(
    df: pd.DataFrame,
    unit: str,
    mt_rank_col: str,
    wt_rank_col: str,
    vaf_col: str,
    tpm_col: Optional[str] = None,
    thr_s1: float = 2.0,
    thr_s3: float = 0.5,
    thr_wt_nonbind: float = 2.0,
) -> Dict[str, Optional[int]]:
    """
    Compute counts for Rank-based S1-S4 under one unit.
    Returns both raw row counts and unique counts (per unit).
    """
    work = df.copy()

    # identifiers
    mt_pep_col = get_col(work, ["MT Epitope Seq", "MT Epitope", "MT_peptide", "MT_PEPTIDE"])
    hla_col = get_col(work, ["HLA Allele", "HLA", "HLA_Allele", "Allele"])

    if mt_pep_col is None:
        raise ValueError("Missing MT peptide column (e.g., 'MT Epitope Seq').")
    if unit == "hla_peptide" and hla_col is None:
        raise ValueError("Missing HLA column (e.g., 'HLA Allele') required for unit=hla_peptide.")

    work["mt_pep"] = work[mt_pep_col].astype(str).fillna("").str.strip()
    work["mutation_id"] = build_mutation_id(work)
    work["hla_norm"] = work[hla_col].apply(normalize_hla) if hla_col is not None else ""

    # numeric
    if mt_rank_col not in work.columns:
        raise ValueError(f"Missing MT rank column: {mt_rank_col}")
    if wt_rank_col not in work.columns:
        raise ValueError(f"Missing WT rank column: {wt_rank_col}")
    if vaf_col not in work.columns:
        raise ValueError(f"Missing VAF column: {vaf_col}")

    work["mt_rank"] = work[mt_rank_col].apply(_safe_float)
    work["wt_rank"] = work[wt_rank_col].apply(_safe_float)
    work["vaf"] = work[vaf_col].apply(_safe_float)

    has_tpm = tpm_col is not None and tpm_col in work.columns
    if has_tpm:
        work["tpm_val"] = work[tpm_col].apply(_safe_float)

    # totals
    total_raw = int(len(work))
    total_unique = int(len(dedup_by_unit(work, unit)))

    # S1
    m1 = work["mt_rank"] < thr_s1

    out: Dict[str, Optional[int]] = {
        "total_raw": total_raw,
        "total_unique": total_unique,
        "s1_rank_raw": int(m1.sum()),
        "s1_rank_unique": int(len(dedup_by_unit(work[m1], unit))),
        "s2_rank_raw": None,
        "s2_rank_unique": None,
        "s3_rank_raw": None,
        "s3_rank_unique": None,
        "s4_rank_raw": None,
        "s4_rank_unique": None,
    }

    if not has_tpm:
        return out

    # S2-S4
    m2 = (work["mt_rank"] < thr_s1) & (work["tpm_val"] > 1)
    m3 = (work["mt_rank"] < thr_s3) & (work["tpm_val"] > 5)
    m4 = m3 & (work["vaf"] > 0.1) & (work["wt_rank"] > thr_wt_nonbind)

    out.update(
        {
            "s2_rank_raw": int(m2.sum()),
            "s2_rank_unique": int(len(dedup_by_unit(work[m2], unit))),
            "s3_rank_raw": int(m3.sum()),
            "s3_rank_unique": int(len(dedup_by_unit(work[m3], unit))),
            "s4_rank_raw": int(m4.sum()),
            "s4_rank_unique": int(len(dedup_by_unit(work[m4], unit))),
        }
    )
    return out


def process_file_rank_long(
    path: str,
    units: List[str],
    mt_rank_col: str,
    wt_rank_col: str,
    vaf_col: str,
    tpm_col: Optional[str],
    thr_s1: float,
    thr_s3: float,
    thr_wt_nonbind: float,
    sample_fullname: bool = False,
) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    base = os.path.basename(path)
    name_no_ext = base.rsplit(".", 1)[0]
    sample = name_no_ext if sample_fullname else name_no_ext.split("_")[0]

    records = []
    for u in units:
        counts = compute_tnb_rank_counts(
            df=df,
            unit=u,
            mt_rank_col=mt_rank_col,
            wt_rank_col=wt_rank_col,
            vaf_col=vaf_col,
            tpm_col=tpm_col,
            thr_s1=thr_s1,
            thr_s3=thr_s3,
            thr_wt_nonbind=thr_wt_nonbind,
        )

        # Convert to long: one row per strategy x metric (raw/unique), plus totals
        # Totals (not strategy-specific)
        records.append(
            {
                "sample": sample,
                "file": base,
                "unit": u,
                "strategy": "total",
                "metric": "raw",
                "value": counts["total_raw"],
            }
        )
        records.append(
            {
                "sample": sample,
                "file": base,
                "unit": u,
                "strategy": "total",
                "metric": "unique",
                "value": counts["total_unique"],
            }
        )

        for strat in ["s1_rank", "s2_rank", "s3_rank", "s4_rank"]:
            raw_k = f"{strat}_raw"
            uniq_k = f"{strat}_unique"
            records.append(
                {
                    "sample": sample,
                    "file": base,
                    "unit": u,
                    "strategy": strat,
                    "metric": "raw",
                    "value": counts.get(raw_k, None),
                }
            )
            records.append(
                {
                    "sample": sample,
                    "file": base,
                    "unit": u,
                    "strategy": strat,
                    "metric": "unique",
                    "value": counts.get(uniq_k, None),
                }
            )

    return pd.DataFrame(records)


def collect_files(input_path: str, pattern: str = "*.tsv") -> List[str]:
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        return sorted(glob.glob(os.path.join(input_path, pattern)))
    raise FileNotFoundError(f"未找到路径: {input_path}")


def main():
    parser = argparse.ArgumentParser(
        description="计算 QC 后文件的 TNB（NetMHCpan %Rank 版本；long-format 输出）。"
    )
    parser.add_argument("input", help="QC 后的 TSV 文件或包含 TSV 的目录")
    parser.add_argument("--pattern", default="*.tsv", help="输入为目录时匹配模式，默认 *.tsv")

    parser.add_argument(
        "--units",
        default="mutation",
        help="计数单位，逗号分隔：mutation,peptide,hla_peptide（默认 mutation）",
    )

    parser.add_argument(
        "--mt-rank-col",
        default="NetMHCpan MT Percentile",
        help="MT %Rank/%Percentile 列名（默认 NetMHCpan MT Percentile）",
    )
    parser.add_argument(
        "--wt-rank-col",
        default="NetMHCpan WT Percentile",
        help="WT %Rank/%Percentile 列名（默认 NetMHCpan WT Percentile）",
    )

    parser.add_argument(
        "--vaf-col",
        default="Tumor DNA VAF",
        help="VAF 列名（默认 Tumor DNA VAF）",
    )
    parser.add_argument(
        "--tpm-col",
        default="TPM",
        help="TPM 列名（默认 TPM；若该列不存在会自动降级只算S1_rank）",
    )

    # thresholds (defaults = what you confirmed)
    parser.add_argument("--thr-s1", type=float, default=2.0, help="S1/S2 的 mt_rank 阈值（默认 2.0）")
    parser.add_argument("--thr-s3", type=float, default=0.5, help="S3/S4 的 mt_rank 阈值（默认 0.5）")
    parser.add_argument(
        "--thr-wt-nonbind",
        type=float,
        default=2.0,
        help="S4 的 wt_rank 不结合阈值（默认 2.0）",
    )

    parser.add_argument(
        "--sample-fullname",
        action="store_true",
        help="使用去扩展名的完整文件名作为sample（默认取下划线前第一段）。",
    )
    parser.add_argument(
        "--output",
        help="输出CSV路径。默认写到输入目录/文件所在目录：TNB_rank_long.csv",
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

    frames = []
    for f in files:
        try:
            frames.append(
                process_file_rank_long(
                    path=f,
                    units=units,
                    mt_rank_col=args.mt_rank_col,
                    wt_rank_col=args.wt_rank_col,
                    vaf_col=args.vaf_col,
                    tpm_col=args.tpm_col if args.tpm_col else None,
                    thr_s1=args.thr_s1,
                    thr_s3=args.thr_s3,
                    thr_wt_nonbind=args.thr_wt_nonbind,
                    sample_fullname=args.sample_fullname,
                )
            )
        except Exception as e:
            print(f"处理文件失败 {f}: {e}")

    if not frames:
        print("没有成功处理任何文件")
        return

    out_df = pd.concat(frames, ignore_index=True)

    # output path
    if args.output:
        out_path = args.output
    else:
        base_dir = os.path.dirname(args.input) if os.path.isfile(args.input) else args.input
        out_path = os.path.join(base_dir, "TNB_rank_long.csv")

    out_df.to_csv(out_path, index=False)
    print(f"已保存 long-format 汇总: {out_path}")
    print(out_df.head(12))


if __name__ == "__main__":
    main()
