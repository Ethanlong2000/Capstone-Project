#!/bin/bash
set -eo pipefail

# ========== 配置 ==========
RAW_DIR="/work/longyh/BY/raw"          # .sralite.1 文件所在目录
FASTQ_OUT_DIR="/work/longyh/BY/fastq/WES"
LOG_DIR="/work/longyh/BY/processed/logs/fasterq"
THREADS_PER_JOB=16
TEMP_DIR="/dev/shm/sra_tmp"

# ========== 初始化 ==========
mkdir -p "$FASTQ_OUT_DIR" "$TEMP_DIR" "$LOG_DIR"

# 依赖检查（确认版本）
if ! command -v fasterq-dump &> /dev/null; then
    echo "[ERROR] fasterq-dump not found" >&2
    exit 1
fi
FASTQ_DUMP_VER=$(fasterq-dump --version | head -1 | awk '{print $2}')
echo "[INFO] Using fasterq-dump version: $FASTQ_DUMP_VER"

# 获取所有 .sralite.1 文件
shopt -s nullglob
sra_files=( "$RAW_DIR"/*.sralite.1 )
shopt -u nullglob

if (( ${#sra_files[@]} == 0 )); then
    echo "[WARN] No .sralite.1 files found in $RAW_DIR"
    exit 0
fi

echo "[INFO] Found ${#sra_files[@]} SRA lite files to process"

# ========== 处理函数（适配 3.2.1 版本） ==========
process_sra() {
    local sra_lite_file="$1"  # 直接传入 .sralite.1 文件完整路径
    local out_dir="$2"
    local threads="$3"
    
    # 提取样本名（去掉 .sralite.1 后缀）
    local base=$(basename "$sra_lite_file" .sralite.1)
    local R1="$out_dir/${base}_1.fastq.gz"
    local R2="$out_dir/${base}_2.fastq.gz"
    local tmpdir="$TEMP_DIR/${base}_$$"
    local log_file="$LOG_DIR/${base}.log"

    # 清空旧日志
    > "$log_file"
    echo "========== [$(date)] START $base ==========" >> "$log_file"

    # 跳过已完成的样本（检查文件存在且非空）
    if [[ -f "$R1" && -f "$R2" && -s "$R1" && -s "$R2" ]]; then
        echo "[SKIP] $base (already completed)" >> "$log_file"
        cat "$log_file"
        return 0
    fi

    # 创建临时目录
    mkdir -p "$tmpdir" || {
        echo "[ERROR] Failed to create temp dir $tmpdir" >> "$log_file"
        echo "========== [$(date)] FAILED $base ==========" >> "$log_file"
        cat "$log_file"
        return 1
    }

    echo "[START] Processing $base (file: $sra_lite_file)" >> "$log_file"

    # 核心修正：适配 fasterq-dump 3.2.1 版本的参数
    # 移除 --resume/--lite/--dir，直接传入 .sralite.1 文件路径
    fasterq-dump \
        --split-3 \
        --skip-technical \
        --threads "$threads" \
        --outdir "$tmpdir" \
        --temp "$tmpdir" \
        --progress \
        "$sra_lite_file" >> "$log_file" 2>&1

    # 检查 fasterq-dump 执行状态
    if [[ $? -ne 0 ]]; then
        echo "[ERROR] fasterq-dump failed for $base (exit code: $?)" >> "$log_file"
        rm -rf "$tmpdir"
        echo "========== [$(date)] FAILED $base ==========" >> "$log_file"
        cat "$log_file"
        return 1
    fi

    # 检查生成的 fastq 文件
    local tmp_R1="$tmpdir/${base}_1.fastq"
    local tmp_R2="$tmpdir/${base}_2.fastq"
    if [[ ! -f "$tmp_R1" || ! -f "$tmp_R2" ]]; then
        echo "[ERROR] Paired files missing: $tmp_R1 or $tmp_R2" >> "$log_file"
        # 列出临时目录内容，方便排查
        ls -l "$tmpdir/" >> "$log_file"
        rm -rf "$tmpdir"
        echo "========== [$(date)] FAILED $base ==========" >> "$log_file"
        cat "$log_file"
        return 1
    fi

    # 压缩为 gz 格式（并行压缩）
    echo "[INFO] Compressing $base R1/R2..." >> "$log_file"
    pigz -p "$threads" -c "$tmp_R1" > "$R1"
    pigz -p "$threads" -c "$tmp_R2" > "$R2"

    # 验证压缩结果
    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "[ERROR] Compression failed (empty output files)" >> "$log_file"
        rm -f "$R1" "$R2"
        rm -rf "$tmpdir"
        echo "========== [$(date)] FAILED $base ==========" >> "$log_file"
        cat "$log_file"
        return 1
    fi

    # 清理临时文件
    rm -rf "$tmpdir"
    echo "[DONE] $base processed successfully" >> "$log_file"
    echo "========== [$(date)] SUCCESS $base ==========" >> "$log_file"
    cat "$log_file"
    return 0
}

# ========== 主循环 ==========
total=${#sra_files[@]}
i=1
for sra in "${sra_files[@]}"; do
    base=$(basename "$sra" .sralite.1)
    echo -e "\n[$i/$total] Starting $base ..."
    echo "Log: $LOG_DIR/${base}.log"

    # 调用处理函数
    if process_sra "$sra" "$FASTQ_OUT_DIR" "$THREADS_PER_JOB"; then
        echo "[$i/$total] [SUCCESS] $base"
    else
        echo "[$i/$total] [FAILURE] $base (check log: $log_file)"
    fi
    ((i++))
done

# 清理空临时目录
find "$TEMP_DIR" -type d -empty -delete

echo -e "\n[SUMMARY] All samples processed."
echo "[OUTPUT] Fastq files: $FASTQ_OUT_DIR"
echo "[LOGS] Log files: $LOG_DIR"