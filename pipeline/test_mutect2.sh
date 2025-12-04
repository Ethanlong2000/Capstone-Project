#!/bin/bash

# ===== 配置变量（根据你的环境修改） =====
SAMPLE_GROUP="Pt8"
TUMOR_BAM="/work/longyh/BY/pipeline/bam/Pt8/tumor.fixed.bam"
NORMAL_BAM="/work/longyh/BY/pipeline/bam/Pt8/normal.fixed.bam"
INTERVAL_FILE="/home/longyh/database/TruSeq_Exome_v1.2_GRCh38_split/0000-scattered.interval_list"
REF_GENOME="/home/longyh/database/GRCh38.d1.vd1.fa"
GERMLINE_RESOURCE="/home/longyh/database/af-only-gnomad.hg38.vcf.gz"
PON="/home/longyh/database/1000g_pon.hg38.vcf.gz"

# 输出路径
OUTPUT_VCF="/work/longyh/BY/pipeline/vcf/Pt8/mutect2/0000-scattered.vcf.gz"
LOG_FILE="/work/longyh/BY/pipeline/logs/Pt8/mutect2/0000-scattered.log"

# Mutect2 参数
THREADS=12
MEM="32G"
TUMOR_ID="${SAMPLE_GROUP}_tumor"
NORMAL_ID="${SAMPLE_GROUP}_normal"

# 创建输出目录
mkdir -p "$(dirname "$OUTPUT_VCF")" "$(dirname "$LOG_FILE")"

# ===== 运行 Mutect2 =====
echo "=== 手动测试 Mutect2（区间: $(basename "$INTERVAL_FILE" .interval_list)）===" > "$LOG_FILE"

gatk --java-options "-Xmx${MEM} -XX:+UseParallelGC" Mutect2 \
    -R "$REF_GENOME" \
    -I "$TUMOR_BAM" -tumor "$TUMOR_ID" \
    -I "$NORMAL_BAM" -normal "$NORMAL_ID" \
    --germline-resource "$GERMLINE_RESOURCE" \
    --panel-of-normals "$PON" \
    -L "$INTERVAL_FILE" \
    -O "$OUTPUT_VCF" \
    --native-pair-hmm-threads "$THREADS" >> "$LOG_FILE" 2>&1

echo "Mutect2 完成。日志见: $LOG_FILE"