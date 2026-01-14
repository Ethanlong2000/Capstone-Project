#!/usr/bin/env bash
set -euo pipefail

# 自动扫描目录下的 .sralite.1 文件，并行 fastq-dump 输出 gzip FASTQ。
# 孤立 reads 输出到 *_3.fastq.gz，保证 _1/_2 配对。
# 自动检测 GNU parallel；没有则用 xargs。
# 为避免联网取 RefSeq/WGS 导致 error code 3，预先关闭 sratoolkit 的远端访问。

SRA_LITE_BASE="/work/longyh/BY/raw"   # 存放 .sralite.1 的目录
OUTDIR="/work/longyh/BY/fastq"        # FASTQ 输出目录
TMP_DIR="/work/longyh/BY/tmp"         # fastq-dump 临时目录
THREADS=8

# 关闭远端仓库（若未安装 vdb-config 或失败则继续）
vdb-config -s /repository/remote/enable=false >/dev/null 2>&1 || true

FASTQ_DUMP=$(command -v fastq-dump || true)
if [[ -z "$FASTQ_DUMP" ]]; then
  echo "[ERROR] 未找到 fastq-dump，请先安装/加载 sratoolkit（例如 conda activate sratoolkit）" >&2
  exit 1
fi

PARALLEL_BIN=$(command -v parallel || true)
PIGZ_BIN=$(command -v pigz || true)

mkdir -p "$OUTDIR" "$TMP_DIR"

# 收集 *.sralite.1 基名作为 RUN
mapfile -d '' RUN_FILES < <(find "$SRA_LITE_BASE" -maxdepth 1 -type f -name "*.sralite.1" -print0)
if [[ ${#RUN_FILES[@]} -eq 0 ]]; then
  echo "[WARN] 在目录 $SRA_LITE_BASE 下未找到 .sralite.1 文件" >&2
  exit 0
fi

RUNS=()
for f in "${RUN_FILES[@]}"; do
  base="$(basename "$f")"
  run="${base%.sralite.1}"
  RUNS+=("$run")
done

TOTAL=${#RUNS[@]}
printf "[INFO] 共需处理 %d 个任务\n" "$TOTAL"

process_one() {
  local RUN="$1"
  local INPUT_PATH="$SRA_LITE_BASE/${RUN}.sralite.1"

  if [[ ! -f "$INPUT_PATH" ]]; then
    echo "[WARN] 未找到输入文件: $INPUT_PATH" >&2
    return 0
  fi

  # 若 fastq.gz 已存在则跳过
  if compgen -G "$OUTDIR/${RUN}"*.fastq.gz > /dev/null; then
    echo "[INFO] 已存在 fastq.gz，跳过 $RUN"
    return 0
  fi

  echo "[INFO] 处理 $RUN，输入: $INPUT_PATH"

  TMPDIR="$TMP_DIR" "$FASTQ_DUMP" \
    --split-3 \
    --gzip \
    --skip-technical \
    -O "$OUTDIR" \
    "$INPUT_PATH"

  # 若为单端输出（fastq-dump 默认名为 RUN.fastq.gz），改名为 RUN_1.fastq.gz
  local single="$OUTDIR/${RUN}.fastq.gz"
  if [[ -f "$single" ]]; then
    local target="$OUTDIR/${RUN}_1.fastq.gz"
    if [[ -e "$target" ]]; then
      echo "[WARN] 目标文件已存在，保留原文件: $target" >&2
    else
      mv "$single" "$target"
      echo "[INFO] 单端输出重命名为: $target"
    fi
  fi
}

export -f process_one
export OUTDIR SRA_LITE_BASE FASTQ_DUMP THREADS

if [[ -n "$PARALLEL_BIN" ]]; then
  echo "[INFO] 检测到 GNU parallel，使用并行 (-j $THREADS)"
  printf '%s\n' "${RUNS[@]}" | "$PARALLEL_BIN" --halt soon,fail=1 -j "$THREADS" process_one {}
else
  echo "[INFO] 未检测到 GNU parallel，使用 xargs (-P $THREADS)"
  printf '%s\n' "${RUNS[@]}" | xargs -n1 -P "$THREADS" bash -c 'process_one "$@"' _
fi

echo "[DONE] 全部处理完成，输出目录: $OUTDIR"