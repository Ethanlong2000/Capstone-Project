#!/usr/bin/env bash
set -euo pipefail

# 目标目录
DST="/work/longyh/BY/processed/naive"

# 本机源目录
LOCAL_SRC="/work/longyh/BY/pipeline/results"

# 远程配置（可用环境变量覆盖）
REMOTE_HOST="${REMOTE_HOST:-192.168.0.108}"
REMOTE_USER="${REMOTE_USER:-longyh}"
# 如远端实际为 /work/longyh/BY/pipeline/result（单数），用 REMOTE_SRC=/work/longyh/BY/pipeline/result 覆盖
REMOTE_SRC="${REMOTE_SRC:-/work/longyh/BY/pipeline/results}"

# 并行度（用于远程拷贝）。默认设为2，避免触发远端 MaxStartups 限制
PARALLEL="${PARALLEL:-2}"

# SSH 复用配置（减少握手）
CONTROL_PATH="${CONTROL_PATH:-$HOME/.ssh/cm-%r@%h:%p}"
SSH_BASE_OPTS=(
  -o BatchMode=yes
  -o ConnectTimeout=5
  -o ServerAliveInterval=30
  -o ServerAliveCountMax=6
  -o ControlMaster=auto
  -o ControlPersist=300
  -o ControlPath="$CONTROL_PATH"
)
SSH_CMD=(ssh "${SSH_BASE_OPTS[@]}")

mkdir -p "$DST"

echo "[INFO] Copying local MHC Class I aggregated files..."
find "$LOCAL_SRC" \
  -type f \
  -path '*/pvacseq_vep/MHC_Class_I/*_tumor.MHC_I.all_epitopes.tsv' \
  -exec cp -t "$DST" {} +

# 远程连通性检查 + 预热控制连接
echo "[INFO] Checking remote host ${REMOTE_USER}@${REMOTE_HOST}..."
if "${SSH_CMD[@]}" "${REMOTE_USER}@${REMOTE_HOST}" 'true' 2>/dev/null; then
  echo "[INFO] Remote reachable. Collecting file list from $REMOTE_SRC ..."

  # 远程查找（仅 MHC Class I），以 null 分隔，避免特殊字符问题
  REMOTE_LIST_CMD="find '$REMOTE_SRC' -type f -path '*/pvacseq_vep/MHC_Class_I/*_tumor.MHC_I.all_epitopes.tsv' -print0"
  # 将远端文件列表落到本地临时文件（null 分隔）
  TMP_LIST="$(mktemp)"
  trap 'rm -f "$TMP_LIST"' EXIT
  "${SSH_CMD[@]}" "${REMOTE_USER}@${REMOTE_HOST}" "$REMOTE_LIST_CMD" > "$TMP_LIST"

  RSYNC_SSH="ssh ${SSH_BASE_OPTS[*]}"

  echo "[INFO] Copying from remote with parallel=${PARALLEL} ..."
  set +e
  xargs -0 -I{} -P "$PARALLEL" rsync -a -e "$RSYNC_SSH" "${REMOTE_USER}@${REMOTE_HOST}:{}" "$DST/" < "$TMP_LIST"
  STATUS=$?
  set -e

  if [[ $STATUS -ne 0 ]]; then
    echo "[WARN] Parallel rsync encountered errors (exit $STATUS). Retrying sequentially..."
    xargs -0 -I{} -P 1 rsync -a -e "$RSYNC_SSH" "${REMOTE_USER}@${REMOTE_HOST}:{}" "$DST/" < "$TMP_LIST"
    echo "[INFO] Sequential retry finished."
  else
    echo "[INFO] Remote copy finished."
  fi
else
  echo "[WARN] Remote host ${REMOTE_HOST} not reachable. Skipping remote copy."
fi

echo "[INFO] Done. Summary:"
echo "Total files in $DST:"
ls -1 "$DST" | wc -l
echo "Sample files:"
ls -1 "$DST" | sort | head -n 10