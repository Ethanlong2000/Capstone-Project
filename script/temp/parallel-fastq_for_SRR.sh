#!/usr/bin/env bash
set -euo pipefail

RUN_LIST="/work/longyh/BY/processed/missing_WES_run_naive.txt"
SRA_BASE="/work/longyh/BY/raw"               # 所有 SRR 文件所在路径（扁平目录）
OUTDIR="/work/longyh/BY/fastq"               # fastq 输出目录
TMPDIR="/work/longyh/BY/tmp"                 # fasterq-dump 临时目录
THREADS=8

mkdir -p "$OUTDIR" "$TMPDIR"

# 读取运行列表（去掉空行和注释行）
mapfile -t RUNS < <(grep -vE '^\s*#' "$RUN_LIST" | sed '/^\s*$/d')
TOTAL=${#RUNS[@]}

if [[ $TOTAL -eq 0 ]]; then
  echo "[WARN] 运行列表为空：$RUN_LIST" >&2
  exit 0
fi

echo "[INFO] 共需处理 $TOTAL 个 SRA 任务"

for idx in "${!RUNS[@]}"; do
  RUN="${RUNS[$idx]}"
  current=$(( idx + 1 ))
  echo "[INFO] ($current/$TOTAL) 正在处理 $RUN"

  SRA_PATH="$SRA_BASE/${RUN}.sra"
  if [[ ! -f "$SRA_PATH" ]]; then
    echo "[WARN] 未找到 SRA 文件: $SRA_PATH" >&2
    continue
  fi

  # 若 fastq.gz 已存在则跳过，避免重复处理
  if compgen -G "$OUTDIR/${RUN}"*.fastq.gz > /dev/null; then
    echo "[INFO] 已存在 fastq.gz，跳过 $RUN"
    continue
  fi

  fasterq-dump \
    --split-files \
    --threads "$THREADS" \
    --temp "$TMPDIR" \
    -O "$OUTDIR" \
    "$SRA_PATH"

  # 压缩输出的 fastq
  pigz -p "$THREADS" "$OUTDIR/${RUN}"*.fastq
done

echo "[DONE] 全部处理完成，输出目录: $OUTDIR"