#!/usr/bin/env bash
set -euo pipefail

# FACETS CNV test script (WES tumor/normal)
# 用途：快速验证 FACETS 在你的数据与靶区上的可用性
# 使用前请填写下面路径（保持绝对路径）

# ======= 必填路径（TODO）=======
TUMOR_BAM="/ABS/PATH/TO/tumor.bam"      # tumor BAM/CRAM（如为CRAM需配置REF_FASTA）
NORMAL_BAM="/ABS/PATH/TO/normal.bam"    # normal BAM/CRAM（如为CRAM需配置REF_FASTA）
TARGET_BED="/ABS/PATH/TO/capture.bed"   # 捕获区BED（WES）
SNP_VCF="/ABS/PATH/TO/dbsnp.vcf.gz"     # 已排序&索引的常见SNP位点（VCF.gz + .tbi）
REF_FASTA="/home/longyh/database/GRCh38.d1.vd1.fa"         # 参考基因组FASTA（如CRAM必须）
OUT_DIR="/ABS/PATH/TO/output"           # 输出目录
SAMPLE_ID="TEST01"                      # 样本ID前缀

# ======= 可选参数（按需调整）=======
THREADS=8
MAPQ=20
BASEQ=20
MAX_DEPTH=10000

# ======= 依赖检查 =======
for bin in snp-pileup Rscript samtools; do
  if ! command -v "$bin" >/dev/null 2>&1; then
    echo "[ERROR] missing dependency: $bin" >&2
    exit 1
  fi
done

mkdir -p "$OUT_DIR"
OUT_PREFIX="$OUT_DIR/$SAMPLE_ID"
PILEUP="$OUT_PREFIX.pileup"
export OUT_DIR SAMPLE_ID PILEUP

# ======= 1) 生成 pileup =======
# 注意：snp-pileup 的选项可按数据质量调整
# 常用参数说明：
#   -q : 读段MAPQ阈值
#   -Q : 碱基质量阈值
#   -P : 读段端点过滤阈值
#   -d : 最大深度
#   -l : BED 限制

snp-pileup \
  -q "$MAPQ" \
  -Q "$BASEQ" \
  -P 20 \
  -r 25,0 \
  -d "$MAX_DEPTH" \
  -l "$TARGET_BED" \
  "$SNP_VCF" \
  "$TUMOR_BAM" \
  "$NORMAL_BAM" \
  "$PILEUP"

# ======= 2) FACETS 主流程 =======
Rscript --vanilla - <<'RSCRIPT'

suppressPackageStartupMessages({
  library(facets)
})

out_dir  <- Sys.getenv("OUT_DIR")
sample_id <- Sys.getenv("SAMPLE_ID")
pileup   <- Sys.getenv("PILEUP")

if (out_dir == "" || sample_id == "" || pileup == "") {
  stop("ENV OUT_DIR/SAMPLE_ID/PILEUP not set")
}

# 读取 pileup
r <- readSnpMatrix(pileup)

# 预处理与分段（参数可按需调整）
xx <- preProcSample(r, cval = 50)
oo <- procSample(xx, cval = 150)

# 拟合 CN 模型
fit <- emcncf(oo)

# 保存结果
saveRDS(list(r = r, xx = xx, oo = oo, fit = fit), file = file.path(out_dir, paste0(sample_id, ".facets.rds")))
write.table(fit$cncf, file = file.path(out_dir, paste0(sample_id, ".facets.cncf.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[OK] FACETS done:\n")
cat("  RDS  :", file.path(out_dir, paste0(sample_id, ".facets.rds")), "\n")
cat("  CNCF :", file.path(out_dir, paste0(sample_id, ".facets.cncf.tsv")), "\n")
RSCRIPT


echo "[DONE] FACETS test completed."
