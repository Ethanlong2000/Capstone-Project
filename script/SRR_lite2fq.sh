#!/bin/bash

# 配置
RAW_DIR="/work/longyh/BY/raw"
FASTQ_OUT_DIR="/work/longyh/BY/fastq/WES"
THREADS=64

# ========== 软件检查 ==========
echo "[INFO] Checking required tools..."
if ! command -v pigz &> /dev/null; then
    echo "[ERROR] pigz not found."
    exit 1
fi

USE_FASTERQ=0
if command -v fasterq-dump &> /dev/null; then
    echo "[INFO] Using fasterq-dump."
    USE_FASTERQ=1
elif command -v fastq-dump &> /dev/null; then
    echo "[INFO] Fasterq-dump unavailable. Using fastq-dump."
    USE_FASTERQ=0
else
    echo "[ERROR] Neither fasterq-dump nor fastq-dump available."
    exit 1
fi

# ========== 准备目录 ==========
mkdir -p "$FASTQ_OUT_DIR"
TMP_DIR="/tmp/sra_fastq_$$"
mkdir -p "$TMP_DIR"

# ========== 获取所有 .1 文件 ==========
shopt -s nullglob
sra_files=( "$RAW_DIR"/*.1 )
shopt -u nullglob

if [[ ${#sra_files[@]} -eq 0 ]]; then
    echo "[WARN] No .1 files in $RAW_DIR"
    exit 0
fi

echo "[INFO] Scanning $RAW_DIR: found ${#sra_files[@]} .1 file(s)."

# ========== 处理每个样本 ==========
processed=0
skipped=0

for sra_file in "${sra_files[@]}"; do
    basename_sra=$(basename "$sra_file" .1)

    # 判断是否已处理（paired 或 single）
    paired_done=0
    single_done=0

    if [[ -f "$FASTQ_OUT_DIR/${basename_sra}_1.fastq.gz" && -f "$FASTQ_OUT_DIR/${basename_sra}_2.fastq.gz" ]]; then
        paired_done=1
    fi

    if [[ -f "$FASTQ_OUT_DIR/${basename_sra}.fastq.gz" ]]; then
        single_done=1
    fi

    if [[ $paired_done -eq 1 || $single_done -eq 1 ]]; then
        echo "[SKIP] $basename_sra already processed. Skipping."
        ((skipped++))
        continue
    fi

    # 执行转换
    ((processed++))
    echo "[PROCESS] ($processed) Converting $basename_sra..."

    if [[ $USE_FASTERQ -eq 1 ]]; then
        fasterq-dump \
            --split-3 \
            --skip-technical \
            --threads "$THREADS" \
            --temp "$TMP_DIR" \
            --outdir "$TMP_DIR" \
            "$sra_file"
    else
        fastq-dump \
            --split-3 \
            --skip-technical \
            --outdir "$TMP_DIR" \
            "$sra_file"
    fi

    # 压缩并移动
    if [[ -f "$TMP_DIR/${basename_sra}_1.fastq" && -f "$TMP_DIR/${basename_sra}_2.fastq" ]]; then
        echo "  -> Paired-end detected. Compressing..."
        pigz -p "$THREADS" -c "$TMP_DIR/${basename_sra}_1.fastq" > "$FASTQ_OUT_DIR/${basename_sra}_1.fastq.gz"
        pigz -p "$THREADS" -c "$TMP_DIR/${basename_sra}_2.fastq" > "$FASTQ_OUT_DIR/${basename_sra}_2.fastq.gz"
    elif [[ -f "$TMP_DIR/${basename_sra}.fastq" ]]; then
        echo "  -> Single-end detected. Compressing..."
        pigz -p "$THREADS" -c "$TMP_DIR/${basename_sra}.fastq" > "$FASTQ_OUT_DIR/${basename_sra}.fastq.gz"
    else
        echo "  [ERROR] Conversion failed for $basename_sra. No FASTQ output."
        continue
    fi

    # 清理临时 FASTQ
    rm -f "$TMP_DIR/${basename_sra}"*.fastq
done

# 清理
rm -rf "$TMP_DIR"

echo "[SUMMARY] Processed: $processed | Skipped (already done): $skipped"
echo "[DONE] FASTQ files stored in: $FASTQ_OUT_DIR"