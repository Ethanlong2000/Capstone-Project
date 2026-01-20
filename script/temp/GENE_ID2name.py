import pandas as pd
import mygene

# 读取 TPM（行名假定为 Entrez Gene ID）
tpm = pd.read_csv('/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm.csv', index_col=0)

scopes = 'entrezgene'
mg = mygene.MyGeneInfo()
res = mg.querymany(
    tpm.index.astype(str).tolist(),
    scopes=scopes,
    fields=['symbol', 'name', 'entrezgene', 'ensembl.gene'],
    species='human',
    as_dataframe=True,
    df_index=True,
)

def pick_symbol(df_row):
    sym = df_row.get('symbol')
    if pd.notna(sym):
        return sym
    name = df_row.get('name')
    return name if pd.notna(name) else 'N/A'

def pick_ensembl_gene(val):
    # 规范化 MyGene 返回的 ensembl.gene 字段为一个 ENSG 编号
    # 可能是 NaN/str/dict/list[dict|str]
    if pd.isna(val):
        return 'N/A'
    try:
        # 字符串：直接返回（若是 ENSG）
        if isinstance(val, str):
            return val if val.startswith('ENSG') else 'N/A'
        # 字典：优先取 'gene' 键
        if isinstance(val, dict):
            g = val.get('gene')
            return g if isinstance(g, str) and g.startswith('ENSG') else 'N/A'
        # 列表/元组：遍历找第一个 ENSG
        if isinstance(val, (list, tuple)):
            for item in val:
                if isinstance(item, str):
                    if item.startswith('ENSG'):
                        return item
                elif isinstance(item, dict):
                    g = item.get('gene')
                    if isinstance(g, str) and g.startswith('ENSG'):
                        return g
            return 'N/A'
    except Exception:
        return 'N/A'
    return 'N/A'

# 处理映射：symbol 和 ensembl_gene_id
if isinstance(res, pd.DataFrame):
    symbol_map = {}
    ensembl_map = {}
    for q, row in res.iterrows():
        q_str = str(q)
        symbol_map[q_str] = pick_symbol(row)
        ensembl_map[q_str] = pick_ensembl_gene(row.get('ensembl.gene'))
else:
    symbol_map = {}
    ensembl_map = {}
    for hit in res:
        q_str = str(hit.get('query'))
        symbol_map[q_str] = hit.get('symbol') or hit.get('name') or 'N/A'
        ens_val = hit.get('ensembl', {}).get('gene') if isinstance(hit.get('ensembl'), dict) else hit.get('ensembl.gene')
        ensembl_map[q_str] = pick_ensembl_gene(ens_val)

# 写入到 TPM 表
tpm['gene_symbol'] = tpm.index.map(lambda x: symbol_map.get(str(x), 'N/A'))
tpm['ensembl_gene_id'] = tpm.index.map(lambda x: ensembl_map.get(str(x), 'N/A'))

# 将新列放到表头
cols = ['gene_symbol', 'ensembl_gene_id'] + [c for c in tpm.columns if c not in ('gene_symbol', 'ensembl_gene_id')]
tpm = tpm[cols]

# 统计打印
total = len(tpm)
symbol_hits = int((tpm['gene_symbol'] != 'N/A').sum())
ensembl_hits = int(tpm['ensembl_gene_id'].str.startswith('ENSG').sum())
symbol_miss = total - symbol_hits
ensembl_miss = total - ensembl_hits

print(f'Total genes: {total}')
print(f'Symbol mapped: {symbol_hits} ({symbol_hits/total:.2%}), missing: {symbol_miss}')
print(f'Ensembl mapped: {ensembl_hits} ({ensembl_hits/total:.2%}), missing: {ensembl_miss}')

# 输出
tpm.to_csv('/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm_with_symbol_ensembl.csv')