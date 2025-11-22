#!/bin/bash

# 定义sratoolkit绝对路径
SRATOOLKIT_BIN="/home/longyh/software/sratoolkit.3.2.1-ubuntu64/bin"
FASTERQ_DUMP="${SRATOOLKIT_BIN}/fasterq-dump"
FASTQ_DUMP="${SRATOOLKIT_BIN}/fastq-dump"

# 定义目录变量
SRA_DIR="/work/longyh/BY/raw/WES"
FASTQ_DIR="/work/longyh/BY/fastq/WES"
TMP_DIR="${FASTQ_DIR}/tmp"  # 临时文件目录
DOWNLOADED_LIST="${SRA_DIR}/downloaded_SRRs.txt"
LOG_FILE="${FASTQ_DIR}/processing.log"
LOCK_FILE="${FASTQ_DIR}/processing.lock"

# 锁机制：防止脚本并行运行
if [[ -f "$LOCK_FILE" ]]; then
    echo "$(date): 已有脚本实例在运行，退出" >> "$LOG_FILE"
    exit 1
fi
touch "$LOCK_FILE"
trap "rm -f '$LOCK_FILE'; rm -rf '${TMP_DIR:?}/'*" EXIT  # 退出时清理锁文件和临时目录

# 时间参数：夜间高资源模式的时间范围（1:00-6:00）
# 强制十进制解析，避免08/09等八进制错误
START_HOUR=$(date +%H -d "1:00")
END_HOUR=$(date +%H -d "6:00")

# 基于服务器配置（80核、755G内存）+ mamba-pigz优化资源参数
MAX_THREADS=32  # 夜间用32核（总核数40%，避免抢占全部资源）
DAY_THREADS=16  # 白天用16核（平衡其他任务）
MAX_MEM="128G"  # 夜间内存限制（总内存17%，避免OOM）
DAY_MEM="64G"   # 白天内存限制
COMPRESS_THREAD_RATIO=5  # 改为整数比例（代表1/2，即×5÷10），避免浮点数

# 初始化目录
mkdir -p "$FASTQ_DIR" "$TMP_DIR"
echo "===== 脚本启动于 $(date) =====" >> "$LOG_FILE"
echo "$(date): 使用sratoolkit路径: ${SRATOOLKIT_BIN}" >> "$LOG_FILE"
echo "$(date): 服务器配置：80核CPU + 755G内存，已优化资源分配" >> "$LOG_FILE"
echo "$(date): 压缩配置：mamba安装的pigz多线程压缩（优先启用）" >> "$LOG_FILE"
echo "$(date): 配置：所有时段处理后均自动压缩（压缩失败不影响处理结果）" >> "$LOG_FILE"
echo "$(date): 已移除--ignore-errors参数（适配3.2.1版本）" >> "$LOG_FILE"

# 检查依赖工具是否存在
if [[ ! -x "$FASTERQ_DUMP" ]]; then
    echo "$(date): fasterq-dump不存在或无执行权限: ${FASTERQ_DUMP}" >> "$LOG_FILE"
    exit 1
fi
if ! command -v pigz &> /dev/null; then
    echo "$(date): 未找到pigz（已安装mamba但未激活环境？），将使用gzip单线程压缩" >> "$LOG_FILE"
    USE_PIGZ=0
else
    echo "$(date): 已检测到pigz（mamba安装），将启用多线程压缩" >> "$LOG_FILE"
    USE_PIGZ=1
fi

# 检查输入文件是否存在
if [[ ! -f "$DOWNLOADED_LIST" ]]; then
    echo "$(date): 未找到下载完成列表 $DOWNLOADED_LIST，脚本退出" >> "$LOG_FILE"
    exit 1
fi

# 批量处理SRA文件
while read -r SRR; do
    [[ -z "$SRR" ]] && continue  # 跳过空行

    # 处理可能带.sra后缀的文件
    SRA_FILE="${SRA_DIR}/${SRR}"
    if [[ ! -f "$SRA_FILE" && -f "${SRA_FILE}.sra" ]]; then
        SRA_FILE="${SRA_FILE}.sra"
    fi

    FASTQ_PREFIX="${FASTQ_DIR}/${SRR}"

    # 检查是否已处理（仅检查压缩格式，确保空间最优）
    if [[ -s "${FASTQ_PREFIX}_1.fastq.gz" || -s "${FASTQ_PREFIX}.fastq.gz" ]]; then
        echo "$(date): 已存在压缩格式文件，跳过: $SRR" >> "$LOG_FILE"
        continue
    fi

    # 检查SRA文件是否有效
    if [[ ! -s "$SRA_FILE" ]]; then
        echo "$(date): SRA文件不存在或为空，跳过: $SRR（文件路径: $SRA_FILE）" >> "$LOG_FILE"
        continue
    fi

    # 获取当前小时（强制十进制）
    current_hour=$(date +%H)
    current_hour_dec=$((10#${current_hour}))
    start_hour_dec=$((10#${START_HOUR}))
    end_hour_dec=$((10#${END_HOUR}))

    # 根据时间分配资源（充分利用硬件，同时留有余地）
    if (( current_hour_dec >= start_hour_dec && current_hour_dec < end_hour_dec )); then
        THREADS=$MAX_THREADS
        MEM_LIMIT=$MAX_MEM
        echo "$(date): 夜间高资源处理: $SRR (线程数: ${THREADS}, 内存限制: ${MEM_LIMIT})" >> "$LOG_FILE"
    else
        THREADS=$DAY_THREADS
        MEM_LIMIT=$DAY_MEM
        echo "$(date): 白天低资源处理: $SRR (线程数: ${THREADS}, 内存限制: ${MEM_LIMIT})" >> "$LOG_FILE"
    fi

    # 修复：整数运算计算压缩线程数（×5÷10 等价于×0.5），并确保至少1线程
    COMPRESS_THREADS=$(( (THREADS * COMPRESS_THREAD_RATIO) / 10 ))
    [[ $COMPRESS_THREADS -lt 1 ]] && COMPRESS_THREADS=1

    # 第一步：尝试用fasterq-dump处理（优先，适配3.2.1版本）
    if "$FASTERQ_DUMP" "$SRA_FILE" \
        -O "$FASTQ_DIR" \
        -t "$TMP_DIR" \
        --split-files \
        --threads "$THREADS" \
        --mem "$MEM_LIMIT" \
        --skip-technical; then
        
        # 标记fasterq-dump处理成功
        echo "$(date): fasterq-dump处理成功: $SRR（开始尝试压缩）" >> "$LOG_FILE"
        
        # 强制压缩（优先mamba-pigz多线程， fallback到gzip）
        # 压缩失败仅记录日志，不删除原文件、不影响流程
        COMPRESS_SUCCESS=0
        if [[ $USE_PIGZ -eq 1 ]]; then
            # mamba-pigz多线程压缩（效率最大化）
            if [[ -f "${FASTQ_PREFIX}_1.fastq" && -f "${FASTQ_PREFIX}_2.fastq" ]]; then
                if pigz -f -p "$COMPRESS_THREADS" "${FASTQ_PREFIX}_1.fastq" "${FASTQ_PREFIX}_2.fastq"; then
                    COMPRESS_SUCCESS=1
                    echo "$(date): 双端文件pigz多线程压缩完成: $SRR（压缩线程数: ${COMPRESS_THREADS}）" >> "$LOG_FILE"
                else
                    echo "$(date): 双端文件pigz压缩失败！保留原fastq文件，可手动重试: $SRR" >> "$LOG_FILE"
                fi
            elif [[ -f "${FASTQ_PREFIX}.fastq" ]]; then
                if pigz -f -p "$COMPRESS_THREADS" "${FASTQ_PREFIX}.fastq"; then
                    COMPRESS_SUCCESS=1
                    echo "$(date): 单端文件pigz多线程压缩完成: $SRR（压缩线程数: ${COMPRESS_THREADS}）" >> "$LOG_FILE"
                else
                    echo "$(date): 单端文件pigz压缩失败！保留原fastq文件，可手动重试: $SRR" >> "$LOG_FILE"
                fi
            else
                echo "$(date): 未找到fasterq-dump生成的fastq文件，跳过压缩: $SRR" >> "$LOG_FILE"
            fi
        else
            # 无pigz，用gzip单线程压缩（兼容模式）
            if [[ -f "${FASTQ_PREFIX}_1.fastq" && -f "${FASTQ_PREFIX}_2.fastq" ]]; then
                if gzip -f "${FASTQ_PREFIX}_1.fastq" "${FASTQ_PREFIX}_2.fastq"; then
                    COMPRESS_SUCCESS=1
                    echo "$(date): 双端文件gzip压缩完成: $SRR" >> "$LOG_FILE"
                else
                    echo "$(date): 双端文件gzip压缩失败！保留原fastq文件，可手动重试: $SRR" >> "$LOG_FILE"
                fi
            elif [[ -f "${FASTQ_PREFIX}.fastq" ]]; then
                if gzip -f "${FASTQ_PREFIX}.fastq"; then
                    COMPRESS_SUCCESS=1
                    echo "$(date): 单端文件gzip压缩完成: $SRR" >> "$LOG_FILE"
                else
                    echo "$(date): 单端文件gzip压缩失败！保留原fastq文件，可手动重试: $SRR" >> "$LOG_FILE"
                fi
            else
                echo "$(date): 未找到fasterq-dump生成的fastq文件，跳过压缩: $SRR" >> "$LOG_FILE"
            fi
        fi

        # 仅记录压缩状态，不影响后续流程
        if [[ $COMPRESS_SUCCESS -eq 1 ]]; then
            echo "$(date): $SRR - 处理+压缩全部完成" >> "$LOG_FILE"
        else
            echo "$(date): $SRR - 处理完成但压缩失败（保留原文件）" >> "$LOG_FILE"
        fi
    
    # 第二步：fasterq-dump失败，降级用fastq-dump（直接输出压缩格式）
    else
        echo "$(date): fasterq-dump处理失败，尝试用fastq-dump降级处理（直接压缩）: $SRR" >> "$LOG_FILE"
        if "$FASTQ_DUMP" "$SRA_FILE" \
            -O "$FASTQ_DIR" \
            --split-files \
            --skip-technical \
            --gzip \
            --readids \
            --threads "$COMPRESS_THREADS"; then  # 使用计算好的压缩线程数
            echo "$(date): fastq-dump直接压缩处理成功: $SRR" >> "$LOG_FILE"
        else
            echo "$(date): fastq-dump也处理失败，放弃: $SRR" >> "$LOG_FILE"
            # 仅清理残留的未压缩文件（避免垃圾文件），不影响主流程
            rm -f "${FASTQ_PREFIX}"*.fastq 2>/dev/null
        fi
    fi

    # 清理当前SRR产生的临时文件（避免累积）
    rm -rf "${TMP_DIR:?}/"*
    rm -rf "${FASTQ_DIR}/${SRR}.tmp" 2>/dev/null

done < "$DOWNLOADED_LIST"

echo "===== 脚本结束于 $(date) =====" >> "$LOG_FILE"