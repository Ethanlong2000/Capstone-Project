#!/bin/bash

FASTQ_DIR="/work/longyh/BY/fastq/WES"
output_dir="/work/longyh/BY/processed/validation_results"
mkdir -p "$output_dir"

LOG="${output_dir}/paired_fq_validation.log"
VALID_LIST="${output_dir}/valid_samples.txt"
INVALID_LIST="${output_dir}/invalid_samples.txt"

> "$LOG"
> "$VALID_LIST"
> "$INVALID_LIST"

echo "[$(date)] 开始验证 FASTQ 配对一致性..." | tee -a "$LOG"

# 获取所有唯一的 SRR 前缀（不含 _1/_2）
for srr in $(ls "$FASTQ_DIR"/*_1.fastq.gz "$FASTQ_DIR"/*_1.fastq 2>/dev/null | sed 's/_1\.fastq.*$//' | sort -u); do
    srr_base=$(basename "$srr")
    r1_gz="${srr}_1.fastq.gz"
    r2_gz="${srr}_2.fastq.gz"
    r1_raw="${srr}_1.fastq"
    r2_raw="${srr}_2.fastq"

    # 优先使用 .gz 文件
    if [[ -f "$r1_gz" && -f "$r2_gz" ]]; then
        use_gz=1
        r1="$r1_gz"
        r2="$r2_gz"
    elif [[ -f "$r1_raw" && -f "$r2_raw" ]]; then
        use_gz=0
        r1="$r1_raw"
        r2="$r2_raw"
    else
        echo "❌ [$srr_base] 缺少 R1 或 R2 文件" | tee -a "$LOG"
        echo "$srr_base" >> "$INVALID_LIST"
        continue
    fi

    # 检查文件是否为空
    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "❌ [$srr_base] R1 或 R2 文件为空" | tee -a "$LOG"
        echo "$srr_base" >> "$INVALID_LIST"
        continue
    fi

    # 提取前 1000 条 reads 的 ID（去除 /1 /2 或 空格后缀）
    if [[ $use_gz -eq 1 ]]; then
        r1_ids=$(zcat "$r1" | head -n 4000 | awk 'NR%4==1 {split($1, a, /[ /]/); print a[1]}')
        r2_ids=$(zcat "$r2" | head -n 4000 | awk 'NR%4==1 {split($1, a, /[ /]/); print a[1]}')
    else
        r1_ids=$(cat "$r1" | head -n 4000 | awk 'NR%4==1 {split($1, a, /[ /]/); print a[1]}')
        r2_ids=$(cat "$r2" | head -n 4000 | awk 'NR%4==1 {split($1, a, /[ /]/); print a[1]}')
    fi

    # 比较 ID 是否完全一致
    if diff -q <(echo "$r1_ids") <(echo "$r2_ids") >/dev/null; then
        echo "✅ [$srr_base] 配对一致" | tee -a "$LOG"
        echo "$srr_base" >> "$VALID_LIST"
    else
        echo "❌ [$srr_base] 配对不一致（read ID 不匹配）" | tee -a "$LOG"
        echo "$srr_base" >> "$INVALID_LIST"

        # 可选：打印前几条不匹配的 ID 用于调试
        echo "    前3条 R1 ID: $(echo "$r1_ids" | head -n 3 | tr '\n' ' ')" >> "$LOG"
        echo "    前3条 R2 ID: $(echo "$r2_ids" | head -n 3 | tr '\n' ' ')" >> "$LOG"
    fi
done

echo "[$(date)] 验证完成！" | tee -a "$LOG"
echo "  有效样本数: $(wc -l < "$VALID_LIST")" | tee -a "$LOG"
echo "  无效样本数: $(wc -l < "$INVALID_LIST")" | tee -a "$LOG"
echo "详细日志: $LOG"
echo "有效样本列表: $VALID_LIST"
echo "需重处理样本列表: $INVALID_LIST"