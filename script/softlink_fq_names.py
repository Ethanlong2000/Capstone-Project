import pandas as pd
import os
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('symlink_creation.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 配置参数
SRR_fq_path = '/work/longyh/BY/fastq/WES'
target_dir = '/work/longyh/BY/pipeline/data'
mapping_dict_file = '/work/longyh/BY/processed/WES_fastq_check_result.xlsx'

def create_soft_link(src, dst):
    """
    创建软链接，处理已存在的链接/文件冲突
    """
    # 确保目标目录存在
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    
    # 如果目标已存在，先删除（区分文件/链接）
    if os.path.exists(dst) or os.path.islink(dst):
        if os.path.islink(dst):
            os.unlink(dst)  # 删除软链接
            logger.warning(f"已删除原有软链接：{dst}")
        else:
            os.remove(dst)  # 删除普通文件（谨慎操作，这里假设是旧文件）
            logger.warning(f"已删除原有文件：{dst}")
    
    try:
        # 创建软链接（使用绝对路径，避免相对路径导致链接失效）
        src_abs = os.path.abspath(src)
        os.symlink(src_abs, dst)
        
        # 验证链接是否有效
        if os.path.exists(dst) and os.readlink(dst) == src_abs:
            logger.info(f"✅ 成功创建软链接：{dst} -> {src_abs}")
            return True
        else:
            logger.error(f"❌ 软链接创建后验证失败：{dst}")
            return False
    except Exception as e:
        logger.error(f"❌ 创建软链接失败 {src} -> {dst}：{str(e)}")
        return False

def main():
    # 读取映射表
    mapping_df = pd.read_excel(mapping_dict_file)
    
    # 重命名sample
    mapping_df['Sample Name'] = mapping_df['Sample Name'].str.replace('_norm', '_normal').str.replace('_pre', '_tumor')
    
    # 去除指定行
    mapping_df = mapping_df[~mapping_df['Sample Name'].isin(['Pt8_tumor', 'Pt8_normal'])]
    
    missing_runs = set()
    success_count = 0
    fail_count = 0
    
    # 遍历每一行
    for index, row in mapping_df.iterrows():
        run_id = str(row['Run']).strip()
        sample_name = str(row['Sample Name']).strip()
        
        # 构建源文件路径
        fq1_path = os.path.join(SRR_fq_path, f"{run_id}_1.fastq.gz")
        fq2_path = os.path.join(SRR_fq_path, f"{run_id}_2.fastq.gz")
        
        # 检查源文件是否存在
        fq1_exists = os.path.exists(fq1_path)
        fq2_exists = os.path.exists(fq2_path)
        
        if fq1_exists and fq2_exists:
            # 构建目标软链接路径
            target_fq1 = os.path.join(target_dir, f"{sample_name}_1.fastq.gz")
            target_fq2 = os.path.join(target_dir, f"{sample_name}_2.fastq.gz")
            
            # 创建两个文件的软链接
            fq1_ok = create_soft_link(fq1_path, target_fq1)
            fq2_ok = create_soft_link(fq2_path, target_fq2)
            
            if fq1_ok and fq2_ok:
                success_count += 1
            else:
                fail_count += 1
                missing_runs.add(run_id)
        else:
            missing_runs.add(run_id)
            missing_files = []
            if not fq1_exists:
                missing_files.append(fq1_path)
            if not fq2_exists:
                missing_files.append(fq2_path)
            logger.warning(f"⚠️ 跳过Run ID {run_id}：缺失源文件 {', '.join(missing_files)}")
            fail_count += 1
    
    # 输出最终统计
    logger.info("="*50)
    logger.info(f"软链接创建完成！成功：{success_count} 对，失败/缺失：{fail_count} 对")
    if missing_runs:
        logger.info(f"失败/缺失的Run ID列表：{', '.join(sorted(missing_runs))}")
    logger.info("详细日志请查看：symlink_creation.log")
    
    # 验证示例（可选）
    logger.info("\n🔍 验证前5个软链接有效性：")
    target_files = [f for f in os.listdir(target_dir) if f.endswith('.fastq.gz')][:5]
    for f in target_files:
        f_path = os.path.join(target_dir, f)
        if os.path.islink(f_path):
            logger.info(f"  {f_path} -> {os.readlink(f_path)} (有效)")
        else:
            logger.info(f"  {f_path} (不是软链接)")

if __name__ == "__main__":
    main()