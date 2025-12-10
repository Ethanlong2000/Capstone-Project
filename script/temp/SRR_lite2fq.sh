#!/bin/bash

# ========== 配置 ==========
RAW_DIR="/work/longyh/BY/raw"
FASTQ_OUT_DIR="/work/longyh/BY/fastq/WES"

# 优化后的资源配置（双端WES适配）
MAX_JOBS=8              # 保持原配置，匹配/dev/shm容量
THREADS_PER_JOB=16      # 保持原线程数，WES数据量大，单任务多线程更高效
TEMP_DIR="/dev/shm/sra_tmp"  # 统一临时目录

# ========== 依赖检查 ==========
for tool in pigz fasterq-dump parallel; do
    if ! command -v "$tool" &> /dev/null; then
        echo "[ERROR] '$tool' not found." >&2
        exit 1
    fi
done

mkdir -p "$FASTQ_OUT_DIR"
mkdir -p "$TEMP_DIR"

# 获取样本：精准匹配 .sralite.1 后缀的文件
shopt -s nullglob
# 只匹配 *.sralite.1 的文件，符合你的真实文件名格式
sra_files=( "$RAW_DIR"/*.sralite.1 )
shopt -u nullglob

if (( ${#sra_files[@]} == 0 )); then
    echo "[WARN] No .sralite.1 files found in $RAW_DIR"
    exit 0
fi

echo "[INFO] Found ${#sra_files[@]} paired-end SRA files (WES data)."
echo "[INFO] Using /dev/shm for temp files (max $MAX_JOBS concurrent, $THREADS_PER_JOB threads each)."
echo "[INFO] Total threads: $((MAX_JOBS * THREADS_PER_JOB)) / 160"

# ========== 处理函数 ==========
process_sra() {
    local sra="$1"
    local out_dir="$2"
    local threads="$3"
    
    # 关键修复：正确提取样本名（去掉 .sralite.1 后缀）
    # 例如：SRR5134828.sralite.1 → SRR5134828
    local base=$(basename "$sra" .sralite.1)
    local R1="$out_dir/${base}_1.fastq.gz"
    local R2="$out_dir/${base}_2.fastq.gz"
    local tmpdir="$TEMP_DIR/$base"

    # 跳过已完成的双端样本（必须同时存在R1和R2）
    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "[SKIP] $base (already completed)"
        rm -rf "$tmpdir"  # 清理残留临时文件
        return 0
    fi

    # 创建样本专属临时目录
    mkdir -p "$tmpdir" || {
        echo "[ERROR] Failed to create temp dir $tmpdir for $base" >&2
        return 1
    }

    echo "[START] Processing $base (source: $sra)"

    # 执行fasterq-dump：保留日志便于排查，双端数据无需调整split参数
    if ! fasterq-dump \
        --split-3 \
        --skip-technical \
        --threads "$threads" \
        --outdir "$tmpdir" \
        --temp "$tmpdir" \
        --progress \
        --resume \  # 断点续传，避免重复耗时
        "$sra"; then
        echo "[ERROR] fasterq-dump failed for $base (check logs above)" >&2
        rm -rf "$tmpdir"
        return 1
    fi

    # 验证双端文件是否生成（WES双端必存在）
    if [[ ! -f "$tmpdir/${base}_1.fastq" || ! -f "$tmpdir/${base}_2.fastq" ]]; then
        echo "[ERROR] Paired-end files missing for $base (R1/R2 not found)" >&2
        rm -rf "$tmpdir"
        return 1
    fi

    # 压缩双端文件（并行压缩，提升效率）
    echo "[INFO] Compressing $base R1/R2..."
    pigz -p "$threads" -c "$tmpdir/${base}_1.fastq" > "$R1" &
    pigz -p "$threads" -c "$tmpdir/${base}_2.fastq" > "$R2" &
    wait  # 等待双端压缩完成

    # 检查压缩结果
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "[ERROR] Compression failed for $base (R1/R2 missing)" >&2
        rm -rf "$tmpdir"
        return 1
    fi

    # 清理临时文件
    rm -rf "$tmpdir"
    echo "[DONE] Processed $base successfully"
}

export -f process_sra
export FASTQ_OUT_DIR TEMP_DIR

# 并行执行：增加失败重试，适配WES大数据量
parallel --bar -j "$MAX_JOBS" --retries 1 process_sra {} "$FASTQ_OUT_DIR" "$THREADS_PER_JOB" ::: "${sra_files[@]}"

# 清理临时目录
rm -rf "$TEMP_DIR"

echo
echo "[SUCCESS] All WES samples processed! Output: $FASTQ_OUT_DIR"