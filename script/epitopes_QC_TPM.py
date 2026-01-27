"""
用于对 pVACseq 结果文件进行 QC，并添加病人特异性 TPM 值

筛选规则：
- 肽段长度：'Peptide Length' ∈ [8, 9, 10, 11]
- 末端突变过滤：'Mutation Position' 可能为多值（如 "1,3"）。若任一突变位置为 1 或等于肽段长度（C 端），则剔除该行。
- 仅保留蛋白编码基因：'Biotype' == 'protein_coding'
- 转录本支持等级：'Transcript Support Level' <= 2
- 选择标准转录本：'Canonical' 为 True

额外清理：
- 删除已知为空的列：['Transcript CDS Flags', 'Tumor RNA Depth', 'Tumor RNA VAF', 'Normal Depth', 'Normal VAF', 'Gene Expression', 'Transcript Expression']
- 删除与研究无关的列：['Index', 'Gene of Interest', 'cterm_7mer_gravy_score', 'max_7mer_gravy_score', 'difficult_n_terminal_residue', 'c_terminal_cysteine', 'c_terminal_proline', 'cysteine_count', 'n_terminal_asparagine', 'asparagine_proline_bond_count']

新增功能：
- 从输入文件名（如 Pt3_tumor.MHC_I.all_epitopes.tsv）提取病人ID（Pt3）
- 从 TPM 表中提取对应病人列和 ensembl_gene_id 列
- 将 Ensembl Gene ID 与 TPM 表匹配，添加 TPM 列（列名 TPM）；若无任何匹配则不添加并文件名不带 TPM 标记；若有匹配则仅保留 TPM 非空行，文件名带 .TPM.QC.tsv

输出：保存到 /work/longyh/BY/processed/QC 目录；若含有效 TPM 数据则文件名加后缀 .TPM.QC.tsv，否则为 .QC.tsv
"""

import os
import sys
import pandas as pd


DEFAULT_INPUT = "/work/longyh/BY/processed/naive/Pt3_tumor.MHC_I.all_epitopes.tsv"
OUTPUT_DIR = "/work/longyh/BY/processed/QC"
TPM_REF = "/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm_with_ensembl_offline.paired.csv"

EMPTY_COLUMNS = [
    "Transcript CDS Flags",
    "Tumor RNA Depth",
    "Tumor RNA VAF",
    "Normal Depth",
    "Normal VAF",
    "Gene Expression",
    "Transcript Expression",
]
IRRELEVANT_COLUMNS = [
    "Index",
    "Gene of Interest",
    "cterm_7mer_gravy_score",
    "max_7mer_gravy_score",
    "difficult_n_terminal_residue",
    "c_terminal_cysteine",
    "c_terminal_proline",
    "cysteine_count",
    "n_terminal_asparagine",
    "asparagine_proline_bond_count",
]


def _safe_int(x):
    """将值安全转换为 int，失败则返回 None。"""
    try:
        if pd.isna(x):
            return None
        return int(str(x).strip())
    except Exception:
        return None


def _safe_float(x):
    """将值安全转换为 float，失败则返回 None。"""
    try:
        if pd.isna(x):
            return None
        return float(str(x).strip())
    except Exception:
        return None


def _parse_positions(cell):
    """
    将 'Mutation Position' 单元格解析为整数列表。
    支持 "1"、"1,3"、" 2, 5 " 等形式。空或不可解析返回空列表。
    """
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    s = str(cell).strip()
    if not s:
        return []
    parts = [p.strip() for p in s.split(",") if p.strip()]
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except Exception:
            pass
    return out


def _is_protein_coding(biotype):
    if pd.isna(biotype):
        return False
    return str(biotype).strip().lower() == "protein_coding"


def _is_canonical(val):
    """将 'Canonical' 字段规范为布尔 True/False。"""
    if pd.isna(val):
        return False
    s = str(val).strip().lower()
    return s in {"true", "yes", "y", "1"}


def _tsl_leq_2(tsl):
    """'Transcript Support Level' 小于等于 2。不可解析视为很大值（不通过）。"""
    n = _safe_int(tsl)
    if n is None:
        return False
    return n <= 2


def _extract_patient_id(path: str) -> str:
    """从输入文件名中提取病人ID，例如 'Pt3_tumor.MHC_I.all_epitopes.tsv' -> 'Pt3'。"""
    base = os.path.basename(path)
    for sep in ["_", "."]:
        if sep in base:
            candidate = base.split(sep)[0].strip()
            if candidate:
                return candidate
    return base.strip()


def _is_terminal_mutation(row) -> bool:
    """
    判断该行是否为末端突变：突变位置包含 1 或包含肽段长度值。
    """
    positions = _parse_positions(row.get("Mutation Position"))
    length = _safe_int(row.get("Peptide Length"))
    if not positions or length is None:
        return False
    for p in positions:
        if p == 1 or p == length:
            return True
    return False



def add_tpm_column(df: pd.DataFrame, patient_id: str) -> tuple[pd.DataFrame, int, bool]:
    """从 TPM 表中提取对应病人列，与 df 左连接，添加 TPM 列。

    返回 (结果DF, 匹配到非空TPM的行数, 是否找到病人列)。
    若未找到病人列或无任何非空TPM，则返回原始 df（不含 TPM 列）。
    """
    tpm_df = pd.read_csv(TPM_REF, dtype=str)

    # 检查 ensembl_gene_id 列是否存在
    if "ensembl_gene_id" not in tpm_df.columns:
        raise KeyError("TPM 文件中缺少 'ensembl_gene_id' 列")

    patient_col_found = patient_id in tpm_df.columns
    if not patient_col_found:
        # 不抛错，后续直接返回原DF
        return df, 0, False

    # 提取所需列
    tpm_subset = tpm_df[["ensembl_gene_id", patient_id]].copy()
    tpm_subset.rename(columns={patient_id: "TPM"}, inplace=True)

    # 确保 Ensembl Gene ID 列存在
    if "Ensembl Gene ID" not in df.columns:
        raise KeyError("输入文件中缺少 'Ensembl Gene ID' 列")

    merged = pd.merge(
        df,
        tpm_subset,
        left_on="Ensembl Gene ID",
        right_on="ensembl_gene_id",
        how="left",
        indicator=False,
    )

    if "ensembl_gene_id" in merged.columns:
        merged = merged.drop(columns=["ensembl_gene_id"])

    merged["TPM"] = merged["TPM"].apply(lambda x: None if pd.isna(x) or str(x).strip() == "" else str(x).strip())
    matched = merged["TPM"].notna().sum()

    if matched == 0:
        if "TPM" in merged.columns:
            merged = merged.drop(columns=["TPM"])
        return merged, 0, True

    merged = merged[merged["TPM"].notna()]
    return merged, int(matched), True


def run_qc(input_path: str, output_dir: str = OUTPUT_DIR) -> str:
    # 读取 TSV（保持为字符串，后续再按需转换）
    df = pd.read_csv(input_path, sep="\t", dtype=str)

    # 删除已知的空列和与研究无关的列
    cols_to_drop = [c for c in EMPTY_COLUMNS + IRRELEVANT_COLUMNS if c in df.columns]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop, errors="ignore")

    original_count = len(df)

    # 肽段长度过滤：8-11
    df["Peptide Length_int"] = df["Peptide Length"].apply(_safe_int)
    df = df[df["Peptide Length_int"].between(8, 11, inclusive="both")]

    # 末端突变过滤：剔除末端（N/C 端）突变条目
    df = df[~df.apply(_is_terminal_mutation, axis=1)]

    # 仅保留蛋白编码基因
    df = df[df["Biotype"].apply(_is_protein_coding)]

    # 转录本支持等级 <= 2
    df = df[df["Transcript Support Level"].apply(_tsl_leq_2)]

    # 选择标准转录本（Canonical 为 True）
    df = df[df["Canonical"].apply(_is_canonical)]

    # 清理临时列
    if "Peptide Length_int" in df.columns:
        df = df.drop(columns=["Peptide Length_int"])

    # 提取病人ID并添加TPM
    patient_id = _extract_patient_id(input_path)
    df, matched_tpm, has_patient_col = add_tpm_column(df, patient_id)

    # 输出目录与文件名（仅在含有非空 TPM 数据时才标记 TPM 字样）
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.basename(input_path)
    stem = base[:-4] if base.lower().endswith(".tsv") else base

    has_tpm_data = matched_tpm > 0
    if not has_tpm_data and "TPM" in df.columns:
        df = df.drop(columns=["TPM"])

    suffix = ".TPM.QC.tsv" if has_tpm_data else ".QC.tsv"
    out_path = os.path.join(output_dir, f"{stem}{suffix}")
    df.to_csv(out_path, sep="\t", index=False)

    print(f"输入: {input_path}")
    print(f"病人ID: {patient_id}")
    print(f"原始行数: {original_count}")
    print(f"过滤后行数: {len(df)}")
    print(f"TPM匹配行数: {matched_tpm}")
    print(f"TPM列有效: {has_tpm_data}")
    if not has_patient_col:
        print(f"TPM 参考文件中未找到病人列: {patient_id}，已跳过 TPM 合并")
    print(f"已保存: {out_path}")
    return out_path


def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT
    if not os.path.exists(input_path):
        print(f"未找到输入文件: {input_path}")
        print("请在命令行提供 TSV 文件路径，例如:\n  python epitopes_QC.py /path/to/pvacseq_result.tsv")
        sys.exit(1)
    try:
        run_qc(input_path)
    except Exception as e:
        print(f"QC 过程中发生错误: {e}")
        sys.exit(2)


if __name__ == "__main__":
    main()