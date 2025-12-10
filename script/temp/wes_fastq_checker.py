import os
import time
import pandas as pd
from datetime import datetime
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('fastq_checker.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)

def check_fastq_files(run_id, fastq_dir):
    """
    检查指定Run对应的fastq.gz文件是否完整存在
    :param run_id: Run编号（如SRR5134763）
    :param fastq_dir: fastq文件存放目录
    :return: 1（都存在）/0（至少一个不存在）
    """
    if pd.isna(run_id) or run_id.strip() == '':
        return 0
    
    # 定义需要检查的两个文件
    file1 = os.path.join(fastq_dir, f"{run_id}_1.fastq.gz")
    file2 = os.path.join(fastq_dir, f"{run_id}_2.fastq.gz")
    
    # 检查文件是否存在且非空
    file1_exists = os.path.exists(file1) and os.path.getsize(file1) > 0
    file2_exists = os.path.exists(file2) and os.path.getsize(file2) > 0
    
    return 1 if (file1_exists and file2_exists) else 0

def process_data():
    """核心处理函数：提取数据、排序、检查文件、生成Excel"""
    # 配置路径
    csv_path = "/work/longyh/BY/processed/WES_meta_filtered.csv"
    fastq_dir = "/work/longyh/BY/fastq/WES"
    output_xlsx = "/work/longyh/BY/processed/WES_fastq_check_result.xlsx"
    
    try:
        # 1. 检查输入文件是否存在
        if not os.path.exists(csv_path):
            logging.error(f"输入CSV文件不存在：{csv_path}")
            return False
        
        if not os.path.exists(fastq_dir):
            logging.error(f"fastq目录不存在：{fastq_dir}")
            return False
        
        # 2. 读取CSV并提取指定列
        logging.info("开始读取CSV文件...")
        df = pd.read_csv(csv_path, encoding='utf-8')
        
        # 检查必要列是否存在
        required_cols = ['Sample Name', 'Run']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logging.error(f"CSV文件缺少必要列：{missing_cols}")
            logging.info(f"可用列名：{df.columns.tolist()}")
            return False
        
        # 3. 提取目标列并去重、去空值
        result_df = df[required_cols].copy()
        result_df = result_df.dropna(subset=required_cols)  # 去除空值行
        result_df = result_df.drop_duplicates()  # 去除重复行
        
        # 4. 按Sample Name排序
        logging.info("按Sample Name排序数据...")
        result_df = result_df.sort_values(by='Sample Name', ascending=True)
        result_df = result_df.reset_index(drop=True)
        
        # 5. 检查fastq文件并添加fq列
        logging.info("开始检查fastq文件...")
        result_df['fq'] = result_df['Run'].apply(lambda x: check_fastq_files(x, fastq_dir))
        
        # 6. 保存为Excel文件
        logging.info(f"保存结果到Excel文件：{output_xlsx}")
        # 处理Excel写入权限问题
        try:
            result_df.to_excel(output_xlsx, index=False, engine='openpyxl')
        except PermissionError:
            # 如果文件被占用，保存为带时间戳的备份文件
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            backup_xlsx = f"/work/longyh/BY/processed/WES_fastq_check_result_{timestamp}.xlsx"
            logging.warning(f"无法写入目标Excel文件（可能被占用），保存到备份文件：{backup_xlsx}")
            result_df.to_excel(backup_xlsx, index=False, engine='openpyxl')
        
        # 输出统计信息
        total_rows = len(result_df)
        fq_exist_count = result_df['fq'].sum()
        logging.info(f"处理完成！总计{total_rows}条记录，其中{fq_exist_count}条记录有完整的fastq文件")
        
        return True
        
    except Exception as e:
        logging.error(f"处理数据时发生错误：{str(e)}", exc_info=True)
        return False

def main_loop():
    """主循环：无限执行，每次执行后休息1小时"""
    logging.info("启动WES fastq文件检查程序（无限循环模式）")
    logging.info("="*50)
    
    while True:
        # 执行一次数据处理
        logging.info("\n开始新一轮数据处理...")
        process_success = process_data()
        
        # 记录本次执行结果
        if process_success:
            logging.info("本轮处理成功完成")
        else:
            logging.error("本轮处理失败")
        
        # 休息1小时（3600秒）
        logging.info(f"\n开始休息1小时，下次执行时间：{datetime.now() + pd.Timedelta(hours=1)}")
        logging.info("如需停止程序，请按 Ctrl+C 中断")
        
        try:
            # 等待3600秒（1小时）
            time.sleep(3600)
        except KeyboardInterrupt:
            # 捕获Ctrl+C中断
            logging.info("\n接收到中断信号，程序即将退出...")
            logging.info("程序正常结束")
            break

if __name__ == "__main__":
    # 检查依赖
    try:
        import openpyxl
    except ImportError:
        logging.error("缺少必要依赖包，请先安装：pip install openpyxl")
        exit(1)
    
    # 启动主循环
    main_loop()