#!/bin/bash
# filepath: /work/longyh/BY/script/temp/validate_sra.sh

set -euo pipefail

SRA_DIR="/work/longyh/BY/raw"
OUT_DIR="/work/longyh/BY/processed/validation_results"
LOG="${OUT_DIR}/sra_validation.log"
OK_LIST="${OUT_DIR}/valid_sra.txt"
BAD_LIST="${OUT_DIR}/invalid_sra.txt"

# 并行任务数（按机器核数调整）
JOBS=4

mkdir -p "$OUT_DIR"
: > "$LOG"
: > "$OK_LIST"
: > "$BAD_LIST"

echo "[$(date)] 开始批量验证 SRA 文件..." | tee -a "$LOG"

export LOG OK_LIST BAD_LIST

validate_one() {
  sra="$1"
  name="$(basename "$sra")"
  echo "[$(date)] 验证 $name ..." | tee -a "$LOG"
  if vdb-validate "$sra" >>"$LOG" 2>&1; then
    echo "✅ $name" | tee -a "$LOG"
    echo "$name" >> "$OK_LIST"
  else
    echo "❌ $name" | tee -a "$LOG"
    echo "$name" >> "$BAD_LIST"
  fi
}

export -f validate_one

find "$SRA_DIR" -maxdepth 1 -name "*.sra" -print0 \
  | xargs -0 -n1 -P "$JOBS" bash -c 'validate_one "$@"' _

echo "[$(date)] 完成。" | tee -a "$LOG"
echo "有效文件数: $(wc -l < "$OK_LIST")" | tee -a "$LOG"
echo "无效文件数: $(wc -l < "$BAD_LIST")" | tee -a "$LOG"
echo "详细日志: $LOG"
echo "有效列表: $OK_LIST"
echo "无效列表: $BAD_LIST"