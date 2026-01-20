import pandas as pd

# 读入 TPM（索引为 Entrez Gene ID）
tpm = pd.read_csv('/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm.csv', index_col=0)
tpm.index = tpm.index.astype(str)

# 1) NCBI gene2ensembl（只保留人类），显式检查列，若缺失则报错提示
g2e_path = '/work/longyh/BY/processed/gene2ensembl.gz'
g2e = pd.read_csv(
    g2e_path,
    sep='\t',
    compression='gzip',
    dtype=str,
    usecols=['#tax_id', 'GeneID', 'Ensembl_gene_identifier'],  # 加快读取、避免无关列
)
missing_cols = [c for c in ['#tax_id', 'GeneID', 'Ensembl_gene_identifier'] if c not in g2e.columns]
if missing_cols:
    raise ValueError(f"Missing columns in gene2ensembl: {missing_cols}. Please re-download {g2e_path}.")

g2e = g2e[g2e['#tax_id'] == '9606'][['GeneID', 'Ensembl_gene_identifier']].dropna()
g2e = g2e.rename(columns={'GeneID': 'entrez_id', 'Ensembl_gene_identifier': 'ensembl_gene_id'})
g2e = g2e.drop_duplicates(subset=['entrez_id'])  # 去重

# 2) HGNC 完整集兜底（包含 Entrez 与 Ensembl 对应）
hgnc_path = '/work/longyh/BY/processed/hgnc_complete_set.txt'
hgnc = pd.read_csv(hgnc_path, sep='\t', dtype=str)

def col(x: str):
    x = x.lower()
    if 'entrez' in x and 'id' in x:
        return 'entrez_id'
    if 'ensembl' in x and 'gene' in x and 'id' in x:
        return 'ensembl_gene_id'
    return None

lower_cols = {c: col(c) for c in hgnc.columns}
need = [c for c in hgnc.columns if lower_cols.get(c) in ('entrez_id', 'ensembl_gene_id')]
if not need or 'entrez_id' not in lower_cols.values() or 'ensembl_gene_id' not in lower_cols.values():
    raise ValueError(f"Cannot find 'entrez_id'/'ensembl_gene_id' columns in HGNC file: {hgnc_path}\nColumns: {hgnc.columns.tolist()}")

hgnc = hgnc[need].copy().rename(columns={c: lower_cols[c] for c in need})
hgnc = hgnc.dropna(subset=['entrez_id', 'ensembl_gene_id']).drop_duplicates(subset=['entrez_id'])

# 3) 先用 gene2ensembl 合并
df = tpm.copy()
df = df.merge(g2e, left_index=True, right_on='entrez_id', how='left')

# 确保合并后有 ensembl_gene_id 列，即便为空也创建
if 'ensembl_gene_id' not in df.columns:
    df['ensembl_gene_id'] = pd.NA

# 4) 对 gene2ensembl 未命中者，用 HGNC 兜底
mask_missing = df['ensembl_gene_id'].isna()
if mask_missing.any():
    fallback = df.loc[mask_missing, ['entrez_id']].merge(hgnc, on='entrez_id', how='left')
    df.loc[mask_missing, 'ensembl_gene_id'] = fallback['ensembl_gene_id'].values

# 5) 统计与列重排
total = len(df)
hits = int(df['ensembl_gene_id'].notna().sum())
miss = total - hits
print(f'Total genes: {total}')
print(f'Ensembl mapped (offline): {hits} ({hits/total:.2%}), missing: {miss}')

df['ensembl_gene_id'] = df['ensembl_gene_id'].fillna('N/A')
cols = ['ensembl_gene_id'] + [c for c in df.columns if c not in ('ensembl_gene_id', 'entrez_id')]
df = df[cols]

df.to_csv('/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm_with_ensembl_offline.csv')