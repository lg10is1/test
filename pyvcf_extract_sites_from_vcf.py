import vcf
import csv
import multiprocessing as mp
import argparse

def load_positions(position_file):
    """读取包含染色体和位置的文件"""
    positions = {}
    with open(position_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            chromosome, position = row
            if chromosome not in positions:
                positions[chromosome] = set()
            positions[chromosome].add(int(position))
    return positions

def filter_vcf_records(chunk, positions):
    """筛选VCF文件中的重合位点"""
    filtered_records = []
    for record in chunk:
        if record.CHROM in positions and record.POS in positions[record.CHROM]:
            filtered_records.append(record)
    return filtered_records

def process_vcf_in_chunks(vcf_file, positions, num_threads, output_vcf):
    """多线程处理VCF文件"""
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)
    
    # 将所有记录分块
    records = list(vcf_reader)
    chunk_size = len(records) // num_threads + 1
    chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]

    # 多线程筛选
    pool = mp.Pool(num_threads)
    results = pool.starmap(filter_vcf_records, [(chunk, positions) for chunk in chunks])

    # 写入VCF文件
    for filtered_records in results:
        for record in filtered_records:
            vcf_writer.write_record(record)

    # 关闭资源
    pool.close()
    pool.join()
    vcf_writer.close()

def main():
    """主程序"""
    # 参数解析
    parser = argparse.ArgumentParser(description="Filter VCF records based on positions.")
    parser.add_argument('-p', '--position_file', required=True, help="Path to the position file.")
    parser.add_argument('-v', '--vcf_file', required=True, help="Path to the input VCF file.")
    parser.add_argument('-o', '--output_vcf', required=True, help="Path to the output VCF file.")
    parser.add_argument('-t', '--threads', type=int, default=4, help="Number of threads to use (default: 4).")
    args = parser.parse_args()

    # 加载位置信息
    positions = load_positions(args.position_file)

    # 多线程处理VCF文件
    process_vcf_in_chunks(args.vcf_file, positions, args.threads, args.output_vcf)
    print(f"Filtered VCF saved to {args.output_vcf}")

if __name__ == "__main__":
    main()