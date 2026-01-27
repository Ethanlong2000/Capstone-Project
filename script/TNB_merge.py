###使用peptide+HLA作为原子事件
# 用于从qc后的文件计算不同策略下的TNB
#文件中的列：['Chromosome', 'Start', 'Stop', 'Reference', 'Variant', 'Transcript',
    #    'Transcript Support Level', 'Transcript Length', 'Canonical',
    #    'MANE Select', 'Biotype', 'Ensembl Gene ID', 'Variant Type', 'Mutation',
    #    'Protein Position', 'Gene Name', 'HGVSc', 'HGVSp', 'HLA Allele',
    #    'Peptide Length', 'Sub-peptide Position', 'Mutation Position',
    #    'MT Epitope Seq', 'WT Epitope Seq', 
    #    'Tumor DNA VAF', 'Median MT IC50 Score', 'Median WT IC50 Score',
    #    'Median Fold Change', 'Median MT Percentile', 'Median WT Percentile',
    #    'NetMHCpan WT IC50 Score', 'NetMHCpan MT IC50 Score']

import argparse
import glob
import os
from typing import Dict, List

import pandas as pd


def _safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def _dedup(df: pd.DataFrame) -> pd.DataFrame:
    """以肽段+HLA为唯一ID去重。"""
    return df.drop_duplicates(subset=["MT Epitope Seq", "HLA Allele"], keep="first")


def compute_tnb_counts(df: pd.DataFrame) -> Dict[str, int]:
    """
    计算四种策略下的 TNB（唯一 peptide+HLA 的条目数量）。
    若缺少 TPM 列，则仅返回策略1，其余为 None。
    """
    # 数值化所需字段（指定使用 NetMHCpan 的 IC50 列、VAF 与 TPM）
    df = df.copy()
    df["mt_ic50"] = df.get("NetMHCpan MT IC50 Score", pd.Series(dtype=float)).apply(_safe_float)
    df["wt_ic50"] = df.get("NetMHCpan WT IC50 Score", pd.Series(dtype=float)).apply(_safe_float)
    df["vaf"] = df.get("Tumor DNA VAF", pd.Series(dtype=float)).apply(_safe_float)

    has_tpm = "TPM" in df.columns
    if has_tpm:
        df["tpm_val"] = df["TPM"].apply(_safe_float)

    # 总体计数（不加任何筛选）
    total_raw = int(len(df))
    try:
        total_unique = int(len(_dedup(df)))
    except Exception:
        # 若关键列缺失，避免报错；此时返回None
        total_unique = None

    # 策略1：Binding-only
    mask1 = df["mt_ic50"] < 500
    result = {
        "total_raw": total_raw,
        "total_unique": total_unique,
        "binding_only_raw": int(mask1.sum()),
        "binding_only_unique": int(len(_dedup(df[mask1]))),
        "ic50_500_tpm1_raw": None,
        "ic50_500_tpm1_unique": None,
        "ic50_50_tpm5_raw": None,
        "ic50_50_tpm5_unique": None,
        "high_quality_raw": None,
        "high_quality_unique": None,
    }

    if not has_tpm:
        return result

    mask2 = (df["mt_ic50"] < 500) & (df["tpm_val"] > 1)
    mask3 = (df["mt_ic50"] < 50) & (df["tpm_val"] > 5)
    mask4 = mask3 & (df["vaf"] > 0.1) & (df["wt_ic50"] > 1000)

    result.update(
        {
            "ic50_500_tpm1_raw": int(mask2.sum()),
            "ic50_500_tpm1_unique": int(len(_dedup(df[mask2]))),
            "ic50_50_tpm5_raw": int(mask3.sum()),
            "ic50_50_tpm5_unique": int(len(_dedup(df[mask3]))),
            "high_quality_raw": int(mask4.sum()),
            "high_quality_unique": int(len(_dedup(df[mask4]))),
        }
    )
    return result


def process_file(path: str) -> Dict[str, int]:
    """读取单个 TSV 文件并计算 TNB 统计。"""
    df = pd.read_csv(path, sep="\t", dtype=str)
    counts = compute_tnb_counts(df)
    base = os.path.basename(path)
    # 样本名仅保留文件名（去扩展名）下划线分割的第一段
    name_no_ext = base.rsplit(".", 1)[0]
    sample = name_no_ext.split("_")[0]
    return {"sample": sample, **counts}


def collect_files(input_path: str, pattern: str = "*.tsv") -> List[str]:
    if os.path.isfile(input_path):
        return [input_path]
    if os.path.isdir(input_path):
        files = sorted(glob.glob(os.path.join(input_path, pattern)))
        return files
    raise FileNotFoundError(f"未找到路径: {input_path}")


def main():
    parser = argparse.ArgumentParser(description="计算 QC 后文件的 TNB 不同策略计数（输出CSV）。")
    parser.add_argument("input", help="QC 后的 TSV 文件或包含 TSV 的目录")
    parser.add_argument(
        "--pattern",
        default="*.tsv",
        help="当输入为目录时，用于匹配文件的 glob 模式，默认 *.tsv",
    )
    parser.add_argument(
        "--output",
        help="汇总结果输出路径（csv）。默认写入输入目录/文件所在目录下的 TNB_summary.csv",
    )
    args = parser.parse_args()

    files = collect_files(args.input, args.pattern)
    if not files:
        print("未找到匹配的文件")
        return

    rows = []
    for f in files:
        try:
            rows.append(process_file(f))
        except Exception as e:
            print(f"处理文件失败 {f}: {e}")

    summary = pd.DataFrame(rows)

    # 默认输出路径（CSV）
    if args.output:
        out_path = args.output
    else:
        base_dir = os.path.dirname(args.input) if os.path.isfile(args.input) else args.input
        out_path = os.path.join(base_dir, "TNB_summary.csv")

    # 总是以CSV格式写出
    summary.to_csv(out_path, index=False)
    print(f"已保存汇总: {out_path}")
    print(summary)


if __name__ == "__main__":
    main()
