#!/bin/bash
set -euo pipefail

# 配置参数
# SRR_LIST="/work/longyh/BY/processed/undownload_SRR_list.txt"
SRR_LIST="/work/longyh/BY/processed/SRR_list.txt"
RAW_DIR="/work/longyh/BY/raw/WES"            # .sra 文件所在目录
FASTQ_DIR="/work/longyh/BY/fastq/WES"
LOG_DIR="/work/longyh/BY/processed/logs/fasterq"
LOG_FILE="$LOG_DIR/parallel_fastq_$(date +%Y%m%d_%H%M%S).log"

THREADS=16          # 根据服务器资源调整
TMP_DIR="/tmp/parallel_fastq_tmp_$(date +%s)_$$"

# 检查命令是否存在
command -v parallel-fastq-dump >/dev/null || { echo "ERROR: parallel-fastq-dump not found" >&2; exit 1; }

# 创建必要目录
mkdir -p "$FASTQ_DIR" "$TMP_DIR" "$LOG_DIR"

# 日志同时输出到终端和文件
exec > >(tee -a "$LOG_FILE") 2>&1

# 打印开始信息
echo "=== Conversion started at $(date) ==="
echo "Using parallel-fastq-dump v$(parallel-fastq-dump -V 2>&1 | head -n1 | awk '{print $NF}')"
echo "Log: $LOG_FILE"
echo "Input dir: $RAW_DIR"
echo "Output dir: $FASTQ_DIR"
echo ""

# 切换到SRA文件目录（核心：让工具能找到.sra文件）
cd "$RAW_DIR" || { echo "ERROR: Failed to enter $RAW_DIR" >&2; exit 1; }

# 遍历SRR列表处理
while IFS= read -r srr; do
    # 跳过空行
    [[ -z "$srr" ]] && continue

    # 提取纯SRR编号（移除可能的.sra后缀）
    srr_base=$(basename "$srr" .sra)
    
    # 检查.sra文件是否存在（当前目录已切换到RAW_DIR）
    if [[ ! -f "${srr_base}.sra" && ! -f "$srr_base" ]]; then
        echo "ERROR: $RAW_DIR/${srr_base}.sra (or $srr_base) not found. Skipping."
        continue
    fi

    echo "[$(date)] Processing: $srr_base"

    # 核心修正：仅使用工具支持的参数，--sra-id指定纯SRR编号
    if parallel-fastq-dump \
        --sra-id "$srr_base" \
        --threads "$THREADS" \
        --outdir "$FASTQ_DIR" \
        --tmpdir "$TMP_DIR" \
        --split-files \
        --gzip; then
        echo "[$(date)] SUCCESS: $srr_base → compressed FASTQ"
    else
        echo "[$(date)] FAILED: $srr_base"
        continue  # 单个失败不终止整体脚本
    fi

done < "$SRR_LIST"

# 切回原目录（可选）
cd - >/dev/null 2>&1

# 清理临时目录
rm -rf "$TMP_DIR"

echo "=== All done at $(date) ==="