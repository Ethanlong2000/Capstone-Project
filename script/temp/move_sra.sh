#!/bin/bash

set -euo pipefail

SOURCE_PATTERN="/home/longyh/SRR*/SRR*.sra"
DEST_DIR="/work/longyh/BY/raw/WES"
LOG="copy_sra_simple.log"

mkdir -p "$DEST_DIR"

echo "[$(date)] 开始复制所有 SRR*.sra 文件（保留 .sra 后缀）..." | tee "$LOG"
echo "源: $SOURCE_PATTERN" | tee -a "$LOG"
echo "目标: $DEST_DIR" | tee -a "$LOG"

# 使用 rsync 直接匹配 shell glob
if rsync -av --partial --progress $SOURCE_PATTERN "$DEST_DIR/" 2>&1 | tee -a "$LOG"; then
    echo "[$(date)] 复制完成！" | tee -a "$LOG"
else
    echo "[$(date)] 复制过程中出现错误（但 rsync 会尽量继续）" | tee -a "$LOG"
fi