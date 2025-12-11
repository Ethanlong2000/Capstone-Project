#!/bin/bash

# 配置路径
FASTQ_DIR="/work/longyh/BY/fastq/WES"
# INVALID_LIST="/work/longyh/BY/processed/validation_results/invalid_samples.txt"
# INVALID_LIST="/work/longyh/BY/processed/SRR_list.txt"
LOG="/work/longyh/BY/processed/validation_results/remove_invalid.log"

# 日志函数
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG"
}

# 检查必要文件是否存在
if [[ ! -f "$INVALID_LIST" ]]; then
    echo "错误: $INVALID_LIST 不存在！" >&2
    exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "错误: FASTQ 目录 $FASTQ_DIR 不存在！" >&2
    exit 1
fi

> "$LOG"  # 清空日志
log "开始删除无效样本的 FASTQ 文件..."

# 逐行读取 invalid_samples.txt
while IFS= read -r srr || [[ -n "$srr" ]]; do
    # 跳过空行和注释行
    [[ -z "$srr" || "$srr" =~ ^[[:space:]]*# ]] && continue

    r1="${FASTQ_DIR}/${srr}_1.fastq.gz"
    r2="${FASTQ_DIR}/${srr}_2.fastq.gz"

    deleted_any=false

    if [[ -f "$r1" ]]; then
        rm -f "$r1"
        log "已删除: $r1"
        deleted_any=true
    else
        log "跳过（不存在）: $r1"
    fi

    if [[ -f "$r2" ]]; then
        rm -f "$r2"
        log "已删除: $r2"
        deleted_any=true
    else
        log "跳过（不存在）: $r2"
    fi

    if [[ "$deleted_any" == "false" ]]; then
        log "⚠️ 警告: 样本 $srr 无任何 FASTQ 文件可删除"
    fi
done < "$INVALID_LIST"

log "删除操作完成。"
echo "详细日志见: $LOG"