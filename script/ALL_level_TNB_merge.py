#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compute TNB under multiple strategies and multiple counting units:
- mutation-level: unique mutation_id
- peptide-level: unique mutation_id + MT peptide (or peptide alone if desired)
- hla_peptide-level: unique MT peptide + HLA allele

Strategies (same as your design):
S1: mt_ic50 < 500
S2: mt_ic50 < 500 and TPM > 1
S3: mt_ic50 < 50 and TPM > 5
S4: mt_ic50 < 50 and TPM > 5 and VAF > 0.1 and wt_ic50 > 1000
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
        x = str(x).strip()
        if x == "" or x.lower() in {"na", "nan", "none"}:
            return float("nan")
        return float(x)
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
    # remove leading "HLA-" if present
    a = a[4:] if a.startswith("HLA-") else a

    # common compact forms like A0201, B0702, C0602
    m = re.match(r"^([ABC])\*?(\d{2,3})(:?)(\d{2,3})$", a)
    if m:
        locus = m.group(1)
        f1 = m.group(2)
        f2 = m.group(4)
        return f"HLA-{locus}*{f1}:{f2}"

    # forms like A*02:01 already
    m2 = re.match(r"^([ABC])\*(\d{2,3}):(\d{2,3})$", a)
    if m2:
        return f"HLA-{m2.group(1)}*{m2.group(2)}:{m2.group(3)}"

    # fallback: return with HLA- prefix if seems like A/B/C locus
    if re.match(r"^[ABC].*", a):
        return "HLA-" + a
    return a


def get_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first existing column name from candidates."""
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

    # fallback
    needed = ["Chromosome", "Start", "Reference", "Variant"]
    for c in needed:
        if c not in df.columns:
            # If even these are missing, return row index-based id (not ideal but prevents crash)
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
    """
    Deduplicate by:
      - mutation: mutation_id
      - peptide: mutation_id + MT Epitope Seq  (recommended; keeps peptides tied to mutation)
      - hla_peptide: MT Epitope Seq + HLA_norm
    """
    if unit == "mutation":
        return df.drop_duplicates(subset=["mutation_id"], keep="first")

    if unit == "peptide":
        # peptide-level usually should be mutation-aware to avoid merging identical peptides from different mutations
        return df.drop_duplicates(subset=["mutation_id", "mt_pep"], keep="first")

    if unit == "hla_peptide":
        return df.drop_duplicates(subset=["mt_pep", "hla_norm"], keep="first")

    raise ValueError(f"Unknown unit: {unit}")


# -----------------------------
# Core computation
# -----------------------------
def compute_tnb_counts(
    df: pd.DataFrame,
    unit: str,
    mt_ic50_col: str,
    wt_ic50_col: str,
    vaf_col: str,
    tpm_col: Optional[str] = None,
) -> Dict[str, Optional[int]]:
    """
    Compute counts for S1-S4 under one unit.
    Returns both raw row counts and unique counts (per unit).
    """

    work = df.copy()

    # core identifiers
    mt_pep_col = get_col(work, ["MT Epitope Seq", "MT Epitope", "MT_peptide", "MT_PEPTIDE"])
    hla_col = get_col(work, ["HLA Allele", "HLA", "HLA_Allele", "Allele"])

    if mt_pep_col is None:
        raise ValueError("Missing MT peptide column (e.g., 'MT Epitope Seq').")
    if hla_col is None:
        # allow peptide/mutation units without HLA, but hla_peptide requires it
        if unit == "hla_peptide":
            raise ValueError("Missing HLA column (e.g., 'HLA Allele').")

    work["mt_pep"] = work[mt_pep_col].astype(str).fillna("").str.strip()
    work["mutation_id"] = build_mutation_id(work)

    if hla_col is not None:
        work["hla_norm"] = work[hla_col].apply(normalize_hla)
    else:
        work["hla_norm"] = ""

    # numeric fields
    work["mt_ic50"] = work[mt_ic50_col].apply(_safe_float) if mt_ic50_col in work.columns else float("nan")
    work["wt_ic50"] = work[wt_ic50_col].apply(_safe_float) if wt_ic50_col in work.columns else float("nan")
    work["vaf"] = work[vaf_col].apply(_safe_float) if vaf_col in work.columns else float("nan")

    has_tpm = tpm_col is not None and tpm_col in work.columns
    if has_tpm:
        work["tpm_val"] = work[tpm_col].apply(_safe_float)

    # total counts
    total_raw = int(len(work))
    total_unique = int(len(dedup_by_unit(work, unit)))

    # masks should exclude NaNs automatically via comparisons -> False on NaN
    mask1 = work["mt_ic50"] < 500

    out = {
        "total_raw": total_raw,
        "total_unique": total_unique,
        "s1_raw": int(mask1.sum()),
        "s1_unique": int(len(dedup_by_unit(work[mask1], unit))),
        "s2_raw": None,
        "s2_unique": None,
        "s3_raw": None,
        "s3_unique": None,
        "s4_raw": None,
        "s4_unique": None,
    }

    if not has_tpm:
        return out

    mask2 = (work["mt_ic50"] < 500) & (work["tpm_val"] > 1)
    mask3 = (work["mt_ic50"] < 50) & (work["tpm_val"] > 5)
    mask4 = mask3 & (work["vaf"] > 0.1) & (work["wt_ic50"] > 1000)

    out.update(
        {
            "s2_raw": int(mask2.sum()),
            "s2_unique": int(len(dedup_by_unit(work[mask2], unit))),
            "s3_raw": int(mask3.sum()),
            "s3_unique": int(len(dedup_by_unit(work[mask3], unit))),
            "s4_raw": int(mask4.sum()),
            "s4_unique": int(len(dedup_by_unit(work[mask4], unit))),
        }
    )
    return out


def process_file(
    path: str,
    units: List[str],
    mt_ic50_col: str,
    wt_ic50_col: str,
    vaf_col: str,
    tpm_col: Optional[str],
    sample_fullname: bool = False,
) -> Dict[str, Optional[int]]:
    df = pd.read_csv(path, sep="\t", dtype=str)

    base = os.path.basename(path)
    name_no_ext = base.rsplit(".", 1)[0]
    # 默认取下划线前第一段（如 Pt10_tumor... -> Pt10）
    sample = name_no_ext if sample_fullname else name_no_ext.split("_")[0]

    row: Dict[str, Optional[int]] = {"sample": sample, "file": base}

    for u in units:
        counts = compute_tnb_counts(
            df=df,
            unit=u,
            mt_ic50_col=mt_ic50_col,
            wt_ic50_col=wt_ic50_col,
            vaf_col=vaf_col,
            tpm_col=tpm_col,
        )
        # prefix columns by unit
        for k, v in counts.items():
            row[f"{u}_{k}"] = v

    # self-check fields (useful for debugging HLA formatting issues)
    if "HLA Allele" in df.columns:
        # how many distinct HLA raw vs normalized (rough proxy for formatting issues)
        raw_unique = df["HLA Allele"].astype(str).fillna("").str.strip().nunique()
        norm_unique = df["HLA Allele"].apply(normalize_hla).nunique()
        row["hla_unique_raw"] = int(raw_unique)
        row["hla_unique_norm"] = int(norm_unique)

    return row


def collect_files(input_path: str, pattern: str = "*.tsv") -> List[str]:
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        return sorted(glob.glob(os.path.join(input_path, pattern)))
    raise FileNotFoundError(f"未找到路径: {input_path}")


def main():
    parser = argparse.ArgumentParser(description="计算 QC 后文件的 TNB（支持多策略+多计数单位，输出CSV）。")
    parser.add_argument("input", help="QC 后的 TSV 文件或包含 TSV 的目录")
    parser.add_argument("--pattern", default="*.tsv", help="输入为目录时匹配模式，默认 *.tsv")
    parser.add_argument(
        "--units",
        default="mutation,peptide,hla_peptide",
        help="计数单位，逗号分隔：mutation,peptide,hla_peptide（默认全算）",
    )
    parser.add_argument(
        "--mt-ic50-col",
        default="NetMHCpan MT IC50 Score",
        help="MT IC50 列名（默认 NetMHCpan MT IC50 Score）",
    )
    parser.add_argument(
        "--wt-ic50-col",
        default="NetMHCpan WT IC50 Score",
        help="WT IC50 列名（默认 NetMHCpan WT IC50 Score）",
    )
    parser.add_argument(
        "--vaf-col",
        default="Tumor DNA VAF",
        help="VAF 列名（默认 Tumor DNA VAF）",
    )
    parser.add_argument(
        "--tpm-col",
        default="TPM",
        help="TPM 列名（默认 TPM；若该列不存在会自动降级只算S1）",
    )
    parser.add_argument(
        "--sample-from-filename",
        action="store_true",
        help="兼容旧参数：从文件名解析sample（已默认启用）。",
    )
    parser.add_argument(
        "--sample-fullname",
        action="store_true",
        help="使用去扩展名的完整文件名作为sample（默认取下划线前第一段）。",
    )
    parser.add_argument(
        "--output",
        help="输出CSV路径。默认写到输入目录/文件所在目录：TNB_summary_by_unit.csv",
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

    rows = []
    for f in files:
        try:
            rows.append(
                process_file(
                    path=f,
                    units=units,
                    mt_ic50_col=args.mt_ic50_col,
                    wt_ic50_col=args.wt_ic50_col,
                    vaf_col=args.vaf_col,
                    tpm_col=args.tpm_col if args.tpm_col else None,
                    sample_fullname=args.sample_fullname,
                )
            )
        except Exception as e:
            print(f"处理文件失败 {f}: {e}")

    summary = pd.DataFrame(rows)

    # output path
    if args.output:
        out_path = args.output
    else:
        base_dir = os.path.dirname(args.input) if os.path.isfile(args.input) else args.input
        out_path = os.path.join(base_dir, "TNB_summary_by_unit.csv")

    summary.to_csv(out_path, index=False)
    print(f"已保存汇总: {out_path}")
    print(summary.head())


if __name__ == "__main__":
    main()
