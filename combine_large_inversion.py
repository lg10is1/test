import os
import glob
import pandas as pd
import argparse
import math

def calculate_reciprocal_overlap(s1, e1, s2, e2):
    """计算两个区间的相互重叠率 (Reciprocal Overlap)"""
    overlap_start, overlap_end = max(s1, s2), min(e1, e2)
    if overlap_start >= overlap_end: return 0
    overlap_len = overlap_end - overlap_start
    len1, len2 = e1 - s1, e2 - s2
    return min(overlap_len / len1, overlap_len / len2)

def parse_gt(gt_str):
    """解析基因型，返回变异等位基因拷贝数(0, 1, 2)"""
    if not gt_str or "./." in gt_str or "." == gt_str: return 0
    gt = gt_str.replace('|', '/').split('/')
    count = 0
    for allele in gt:
        if allele == '1': count += 1
    return count

def load_sample_list(file_path):
    """从文本文件加载样本名列表"""
    if not os.path.exists(file_path): return set()
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def hypergeom_pmf(k, M, n, N):
    """计算超几何分布的概率质量函数 (PMF)"""
    try:
        return (math.comb(n, k) * math.comb(M - n, N - k)) / math.comb(M, N)
    except (ValueError, OverflowError):
        return 0.0

def fisher_exact_onesided_greater(a, b, c, d):
    """手动实现 Fisher 精确检验 (单侧富集)"""
    N = a + b + c + d
    n_cases = a + b
    n_carriers = a + c
    p_val = 0.0
    for i in range(a, min(n_cases, n_carriers) + 1):
        p_val += hypergeom_pmf(i, N, n_carriers, n_cases)
    return min(p_val, 1.0)

def calculate_fdr(p_values):
    """手动实现 BH 法 FDR 矫正"""
    n = len(p_values)
    if n == 0: return []
    indexed_p = sorted(enumerate(p_values), key=lambda x: x[1])
    fdr = [0] * n
    for i, (idx, p) in enumerate(indexed_p):
        fdr[idx] = p * n / (i + 1)
    for i in range(n - 2, -1, -1):
        fdr[indexed_p[i][0]] = min(fdr[indexed_p[i][0]], fdr[indexed_p[i+1][0]])
    return [min(f, 1.0) for f in fdr]

def main():
    parser = argparse.ArgumentParser(description="Large Inversion 合并、统计与 Singleton 过滤工具")
    parser.add_argument("-i", "--input_dir", required=True, help="VCF 目录")
    parser.add_argument("-case", "--case_file", required=True, help="case 列表")
    parser.add_argument("-ctrl", "--ctrl_file", required=True, help="control 列表")
    parser.add_argument("-d", "--dist", type=int, default=1000, help="断点距离阈值")
    parser.add_argument("-r", "--overlap", type=float, default=0.5, help="重叠率阈值")
    parser.add_argument("-o", "--output", default="INV_Stats_Full.txt", help="主输出文件名")
    args = parser.parse_args()

    # 1. 加载样本
    case_list = load_sample_list(args.case_file)
    ctrl_list = load_sample_list(args.ctrl_file)
    case_set, ctrl_set = set(case_list), set(ctrl_list)

    # 2. 读取 VCF
    vcf_files = sorted(glob.glob(os.path.join(args.input_dir, "*_filtered.vcf")))
    all_variants = []
    found_samples = []

    print(f"[*] 正在读取 VCF 数据...")
    for vcf_path in vcf_files:
        sample_name = os.path.basename(vcf_path).replace("_filtered.vcf", "")
        if sample_name not in case_set and sample_name not in ctrl_set: continue
        found_samples.append(sample_name)
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith("#"): continue
                cols = line.strip().split('\t')
                pos = int(cols[1])
                info, gt = cols[7], cols[9].split(':')[0]
                end = pos; svlen = 0
                for item in info.split(';'):
                    if item.startswith("END="): end = int(item.split('=')[1])
                    if item.startswith("SVLEN="): svlen = abs(int(item.split('=')[1]))
                if end == pos: end = pos + svlen
                all_variants.append({'chrom': cols[0], 'start': pos, 'end': end, 'sample': sample_name, 'gt': gt})

    # 3. 聚类合并
    all_variants.sort(key=lambda x: (x['chrom'], x['start']))
    merged_clusters = []
    print(f"[*] 正在执行变异聚类 (Dist <= {args.dist}, Overlap >= {args.overlap})...")
    for var in all_variants:
        match_found = False
        for cluster in reversed(merged_clusters[-100:]):
            if var['chrom'] == cluster['chrom'] and \
               abs(var['start'] - cluster['start']) <= args.dist and \
               abs(var['end'] - cluster['end']) <= args.dist and \
               calculate_reciprocal_overlap(var['start'], var['end'], cluster['start'], cluster['end']) >= args.overlap:
                cluster['samples'][var['sample']] = var['gt']
                match_found = True
                break
        if not match_found:
            merged_clusters.append({'chrom': var['chrom'], 'start': var['start'], 'end': var['end'], 'samples': {var['sample']: var['gt']}})

    # 4. 统计与分析
    print(f"[*] 正在计算组间差异及过滤单样本位点...")
    final_data = []
    actual_case = [s for s in case_list if s in found_samples]
    actual_ctrl = [s for s in ctrl_list if s in found_samples]
    N_CASE, N_CTRL = len(actual_case), len(actual_ctrl)

    for c in merged_clusters:
        # 提取各组携带者
        case_gt_list = [parse_gt(c['samples'].get(s, "0/0")) for s in actual_case]
        ctrl_gt_list = [parse_gt(c['samples'].get(s, "0/0")) for s in actual_ctrl]
        
        case_carriers = sum(1 for g in case_gt_list if g > 0)
        ctrl_carriers = sum(1 for g in ctrl_gt_list if g > 0)
        total_carriers = case_carriers + ctrl_carriers  # 用于判定是否为 Singleton

        case_ac = sum(case_gt_list)
        ctrl_ac = sum(ctrl_gt_list)
        
        # Fisher 检验与 OR
        p_val = fisher_exact_onesided_greater(case_carriers, N_CASE - case_carriers, ctrl_carriers, N_CTRL - ctrl_carriers)
        odds_ratio = ((case_carriers + 0.5) * (N_CTRL - ctrl_carriers + 0.5)) / ((N_CASE - case_carriers + 0.5) * (ctrl_carriers + 0.5))

        row = {
            'CHROM': c['chrom'], 'START': c['start'], 'END': c['end'], 'SVLEN': c['end']-c['start'],
            'Total_Carriers': total_carriers,
            'CASE_AC': case_ac, 'CASE_AF': round(case_ac/(N_CASE*2), 6) if N_CASE>0 else 0,
            'CTRL_AC': ctrl_ac, 'CTRL_AF': round(ctrl_ac/(N_CTRL*2), 6) if N_CTRL>0 else 0,
            'OR': round(odds_ratio, 4), 'P_value': p_val
        }
        for s in (actual_case + actual_ctrl): row[s] = c['samples'].get(s, "0/0")
        final_data.append(row)

    # 5. 生成结果与过滤
    df = pd.DataFrame(final_data)
    if df.empty: return
    
    df['FDR_P'] = calculate_fdr(df['P_value'].tolist())
    
    # 调整列顺序
    cols = list(df.columns)
    p_idx = cols.index('P_value')
    df = df[cols[:p_idx+1] + ['FDR_P'] + cols[p_idx+1:-1]]

    # 保存全量结果
    df.sort_values(['P_value', 'CHROM', 'START']).to_csv(args.output, sep='\t', index=False)
    
    # 筛选并保存非单体位点 (Carrier > 1)
    df_no_singletons = df[df['Total_Carriers'] > 1]
    no_singleton_out = args.output.replace(".txt", "_no_singletons.txt")
    df_no_singletons.sort_values(['P_value', 'CHROM', 'START']).to_csv(no_singleton_out, sep='\t', index=False)

    print(f"[#] 分析完成！")
    print(f"    - 原始位点总数: {len(df)}")
    print(f"    - 过滤 Singleton 后位点数 (Carrier > 1): {len(df_no_singletons)}")
    print(f"    - 全量结果: {args.output}")
    print(f"    - 过滤结果: {no_singleton_out}")

if __name__ == "__main__":
    main()