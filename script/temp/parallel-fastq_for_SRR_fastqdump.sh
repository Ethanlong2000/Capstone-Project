#!/usr/bin/env bash
set -euo pipefail

# 使用 fastq-dump 生成 FASTQ（含 gzip 压缩），并将孤立 reads 放在 *_3.fastq.gz，保证 _1/_2 严格配对。
# 自动检测 GNU parallel；没有则用 xargs 并行/串行。检测不到 fastq-dump 时直接报错退出。

RUN_LIST="/work/longyh/BY/processed/missing_WES_run_naive.txt"
SRA_BASE="/work/longyh/BY/raw"    # SRA 扁平目录
OUTDIR="/work/longyh/BY/fastq"    # FASTQ 输出目录
THREADS=8

FASTQ_DUMP=$(command -v fastq-dump || true)
if [[ -z "$FASTQ_DUMP" ]]; then
  echo "[ERROR] 未找到 fastq-dump，请先安装/加载 sratoolkit（例如 conda activate sratoolkit）" >&2
  exit 1
fi

PARALLEL_BIN=$(command -v parallel || true)
PIGZ_BIN=$(command -v pigz || true)

mkdir -p "$OUTDIR"

# 读取运行列表（去掉空行和注释行）
mapfile -t RUNS < <(grep -vE '^\s*#' "$RUN_LIST" | sed '/^\s*$/d')
TOTAL=${#RUNS[@]}

if [[ $TOTAL -eq 0 ]]; then
  echo "[WARN] 运行列表为空：$RUN_LIST" >&2
  exit 0
fi

echo "[INFO] 共需处理 $TOTAL 个 SRA 任务"

process_one() {
  local RUN="$1"
  local SRA_PATH="$SRA_BASE/${RUN}.sra"

  if [[ ! -f "$SRA_PATH" ]]; then
    echo "[WARN] 未找到 SRA 文件: $SRA_PATH" >&2
    return 0
  fi

  # 若 fastq.gz 已存在则跳过
  if compgen -G "$OUTDIR/${RUN}"*.fastq.gz > /dev/null; then
    echo "[INFO] 已存在 fastq.gz，跳过 $RUN"
    return 0
  fi

  echo "[INFO] 处理 $RUN"

  # --split-3: 将未成对 reads 输出到 *_3.fastq.gz，避免扰乱 _1/_2 配对
  "$FASTQ_DUMP" \
    --split-3 \
    --gzip \
    --skip-technical \
    -O "$OUTDIR" \
    "$SRA_PATH"

  # 如需改用 pigz 手动压缩（先去掉 --gzip），可启用下方代码：
  # if [[ -n "$PIGZ_BIN" ]]; then
  #   find "$OUTDIR" -maxdepth 1 -type f -name "${RUN}*.fastq" -print0 | xargs -0 -P1 "$PIGZ_BIN" -p "$THREADS"
  # fi
}

export -f process_one
export OUTDIR SRA_BASE FASTQ_DUMP THREADS

if [[ -n "$PARALLEL_BIN" ]]; then
  echo "[INFO] 检测到 GNU parallel，使用并行 (-j $THREADS)"
  printf '%s\n' "${RUNS[@]}" | "$PARALLEL_BIN" --halt soon,fail=1 -j "$THREADS" process_one {}
else
  echo "[INFO] 未检测到 GNU parallel，使用 xargs (-P $THREADS)"
  printf '%s\n' "${RUNS[@]}" | xargs -n1 -P "$THREADS" bash -c 'process_one "$@"' _
fi

echo "[DONE] 全部处理完成，输出目录: $OUTDIR"
