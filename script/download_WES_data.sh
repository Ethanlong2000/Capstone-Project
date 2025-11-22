#!/bin/bash

SRR_LIST="/work/longyh/BY/processed/WES_download_run_list.txt"
DEST_DIR="/work/longyh/BY/raw/WES"
DOWNLOADED_LIST="${DEST_DIR}/downloaded_SRRs.txt"

mkdir -p "$DEST_DIR"
touch "$DOWNLOADED_LIST"

while read -r SRR; do
    [[ -z "$SRR" ]] && continue

    URL="https://sra-pub-run-odp.s3.amazonaws.com/sra/${SRR}/${SRR}"
    OUTFILE="${DEST_DIR}/${SRR}"

    # 优先检查记录文件
    if grep -q "^${SRR}$" "$DOWNLOADED_LIST"; then
        echo "已记录完成，跳过: $SRR"
        continue
    fi

    # 检查文件是否存在且完整（通过远程大小校验）
    if [[ -f "$OUTFILE" ]]; then
        # 获取远程文件大小（处理代理和可能的换行符）
        remote_size=$(proxychains4 wget --spider --server-response "$URL" 2>&1 | grep -i "Content-Length" | awk '{print $2}' | tr -d '\r')
        # 获取本地文件大小
        local_size=$(stat -c %s "$OUTFILE")

        # 校验大小是否一致
        if [[ -n "$remote_size" && "$local_size" -eq "$remote_size" ]]; then
            echo "文件存在且完整但未记录，补充记录: $SRR"
            echo "$SRR" >> "$DOWNLOADED_LIST"
            continue
        else
            echo "文件存在但不完整（远程大小: $remote_size，本地大小: $local_size），将重新下载: $SRR"
            rm -f "$OUTFILE"  # 删除不完整文件，避免续传失败
        fi
    fi

    # 执行下载
    echo "正在下载: $SRR"
    proxychains4 wget -c --limit-rate=10m -O "$OUTFILE" "$URL"
    if [[ $? -eq 0 ]]; then
        # 下载完成后再次校验大小（双重保险）
        remote_size=$(proxychains4 wget --spider --server-response "$URL" 2>&1 | grep -i "Content-Length" | awk '{print $2}' | tr -d '\r')
        local_size=$(stat -c %s "$OUTFILE")
        if [[ -n "$remote_size" && "$local_size" -eq "$remote_size" ]]; then
            echo "下载完成且校验通过: $SRR，暂停3分钟..."
            echo "$SRR" >> "$DOWNLOADED_LIST"
            sleep 180
        else
            echo "下载完成但校验失败（大小不匹配），未记录: $SRR"
            rm -f "$OUTFILE"  # 删除不完整文件，下次重新下载
        fi
    else
        echo "下载失败: $SRR"
    fi
done < "$SRR_LIST"