set -euo pipefail

GROUP="Pt8"  # 替换为你的样本组名
MEM="32G"

# 路径检查
test -f "bam/${GROUP}/tumor.sorted.dedup.rg.bam"
test -f "bam/${GROUP}/normal.sorted.dedup.rg.bam"

# 修改肿瘤 BAM 的 RG（关键：RGSM 不同于正常）
gatk --java-options "-Xmx${MEM} -XX:+UseParallelGC" AddOrReplaceReadGroups \
  -I "bam/${GROUP}/tumor.sorted.dedup.rg.bam" \
  -O "bam/${GROUP}/tumor.fixed.bam" \
  -RGID "${GROUP}_tumor_RG001" \
  -RGSM "${GROUP}_tumor" \
  -RGLB "${GROUP}_tumor_Lib" \
  -RGPL "ILLUMINA" \
  -RGPU "${GROUP}_tumor_PU001"

samtools index "bam/${GROUP}/tumor.fixed.bam"

# 修改正常 BAM 的 RG
gatk --java-options "-Xmx${MEM} -XX:+UseParallelGC" AddOrReplaceReadGroups \
  -I "bam/${GROUP}/normal.sorted.dedup.rg.bam" \
  -O "bam/${GROUP}/normal.fixed.bam" \
  -RGID "${GROUP}_normal_RG001" \
  -RGSM "${GROUP}_normal" \
  -RGLB "${GROUP}_normal_Lib" \
  -RGPL "ILLUMINA" \
  -RGPU "${GROUP}_normal_PU001"

samtools index "bam/${GROUP}/normal.fixed.bam"

echo "=== 验证肿瘤 BAM 的 SM 信息 ==="
samtools view -H "bam/${GROUP}/tumor.fixed.bam" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/){split($i,a,":"); print a[2]}}'

echo -e "\n=== 验证正常 BAM 的 SM 信息 ==="
samtools view -H "bam/${GROUP}/normal.fixed.bam" | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/){split($i,a,":"); print a[2]}}'