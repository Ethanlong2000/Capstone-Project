import pandas as pd
import mygene

# 读取 TPM
tpm = pd.read_csv('/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm.csv', index_col=0)

# === 可选：小批量抽查行名是否为 Entrez ID（默认注释，需检查可取消注释） ===
# mg_check = mygene.MyGeneInfo()
# sample_ids = tpm.index.astype(str).tolist()[:5]
# print(mg_check.querymany(sample_ids, scopes='entrezgene', fields='symbol,name', species='human'))

# 确认是 Entrez 后，使用单一 scopes
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

# 处理映射：优先 symbol，缺失则 N/A
if isinstance(res, pd.DataFrame):
    symbol_map = res['symbol'].fillna(res.get('name')).to_dict()
else:
    symbol_map = {str(hit['query']): hit.get('symbol') or 'N/A' for hit in res}

tpm['gene_symbol'] = tpm.index.map(lambda x: symbol_map.get(str(x), 'N/A'))

# 将 gene_symbol 放到首列
cols = ['gene_symbol'] + [c for c in tpm.columns if c != 'gene_symbol']
tpm = tpm[cols]

# 输出
tpm.to_csv('/work/longyh/BY/processed/GSE91061_BMS038109Sample.hg19KnownGene.tpm_with_symbol.csv')