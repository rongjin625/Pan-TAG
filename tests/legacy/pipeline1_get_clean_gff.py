import sys

def prepare_gff_filtered(input_gff, output_gff, map_out):
    genes = []
    with open(input_gff, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            if parts[2] != "gene": continue

            attr = parts[8]
            # 过滤逻辑：蛋白编码 且 非通读
            if "gene_type=protein_coding" in attr and "readthrough_gene" not in attr:
                info = {}
                for a in attr.split(';'):
                    if '=' in a:
                        k, v = a.split('=', 1)
                        info[k] = v

                gid = info.get('ID') or info.get('gene_id')
                gname = info.get('gene_name') or gid

                if gid:
                    genes.append({
                        'chrom': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'gid': gid,
                        'gname': gname
                    })

    # 1. 核心排序：染色体升序 -> 起始位置升序 -> 终止位置降序
    # 这样对于同一位置，范围最大的基因（容器）会排在最前面
    genes.sort(key=lambda x: (x['chrom'], x['start'], -x['end']))

    # 2. 嵌套过滤（套娃处理）
    filtered_genes = []
    if genes:
        # 初始放入第一个基因
        last_gene = genes[0]
        filtered_genes.append(last_gene)

        for i in range(1, len(genes)):
            curr = genes[i]

            # 如果是新染色体，直接通过
            if curr['chrom'] != last_gene['chrom']:
                filtered_genes.append(curr)
                last_gene = curr
                continue

            # 由于已经排序：curr['start'] >= last_gene['start'] 是一定的
            # 如果 curr['end'] <= last_gene['end']，说明 curr 完全被包裹在 last_gene 内部
            if curr['end'] <= last_gene['end']:
                # 发现套娃！打印一条信息到 stderr 方便你检查
                sys.stderr.write(f"[Nested] 跳过被包含基因: {curr['gid']} (在 {last_gene['gid']} 内部)\n")
                continue
            else:
                # 只有当 current 的终止点超过了 last 的终止点，它才是一个“新的”位点
                filtered_genes.append(curr)
                last_gene = curr

    # 3. 输出结果
    with open(output_gff, 'w') as f_gff, open(map_out, 'w') as f_map:
        f_map.write("Gene_ID\tGene_Name\n")
        for g in filtered_genes:
            f_gff.write(f"{g['chrom']}\t{g['gid']}\t{g['start']}\t{g['end']}\n")
            f_map.write(f"{g['gid']}\t{g['gname']}\n")

if __name__ == "__main__":
    prepare_gff_filtered(sys.argv[1], sys.argv[2], sys.argv[3])