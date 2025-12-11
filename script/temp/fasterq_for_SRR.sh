#!/bin/bash

set -euo pipefail

SRR_LIST="/work/longyh/BY/processed/SRR_list.txt"
RAW_DIR="/work/longyh/BY/raw/WES"
FASTQ_DIR="/work/longyh/BY/fastq/WES"

# fasterq-dump 线程与内存（单位 MB），建议根据机器内存调整；64G → 64000MB
FASTERQ_THREADS=32
FASTERQ_MEM_MB=64000

# fasterq-dump 使用的临时目录（避免在高占用分区写临时文件）
TMP_DIR="/tmp/fasterq_tmp"

# 压缩线程数（建议 4–8，pigz 收益在 4–8 线程后趋于平缓）
PIGZ_THREADS=8

mkdir -p "$FASTQ_DIR" "$TMP_DIR"

# 检查 pigz 是否可用
if ! command -v pigz &> /dev/null; then
    echo "ERROR: pigz not found. Please install it first:"
    echo "  conda install -c conda-forge pigz"
    echo "  or"
    echo "  sudo apt-get install pigz"
    exit 1
fi

while IFS= read -r srr; do
    [[ -z "$srr" ]] && continue

    SRA_FILE="$RAW_DIR/$srr"
    if [[ ! -f "$SRA_FILE" ]]; then
        echo "ERROR: $SRA_FILE not found! Skipping $srr."
        continue
    fi

    echo "[$(date)] Starting conversion for $srr ..."

    # Step 1: fasterq-dump 转换（不压缩）
    if ! fasterq-dump \
        --split-3 \
        --skip-technical \
        --threads "$FASTERQ_THREADS" \
        --mem "$FASTERQ_MEM_MB" \
        --temp "$TMP_DIR" \
        --outdir "$FASTQ_DIR" \
        "$SRA_FILE"; then
        echo "ERROR: fasterq-dump failed for $srr (see above). Skipping compression."
        continue
    fi

    # Step 2: 用 pigz 并行压缩
    if [[ -f "$FASTQ_DIR/${srr}_1.fastq" ]]; then
        # 双端测序
        echo "  -> Compressing paired-end files with pigz (-p $PIGZ_THREADS)"
        pigz -f -p "$PIGZ_THREADS" "$FASTQ_DIR/${srr}_1.fastq"
        pigz -f -p "$PIGZ_THREADS" "$FASTQ_DIR/${srr}_2.fastq"
    elif [[ -f "$FASTQ_DIR/${srr}.fastq" ]]; then
        # 单端测序
        echo "  -> Compressing single-end file with pigz (-p $PIGZ_THREADS)"
        pigz -f -p "$PIGZ_THREADS" "$FASTQ_DIR/${srr}.fastq"
    else
        echo "WARNING: No FASTQ output found for $srr after fasterq-dump."
    fi

    echo "[$(date)] Completed $srr"
done < "$SRR_LIST"

echo "[$(date)] All samples processed and compressed with pigz."