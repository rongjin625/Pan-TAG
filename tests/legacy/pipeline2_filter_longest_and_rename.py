import sys
from Bio import SeqIO

def rename_and_filter_strict(raw_faa, filtered_gff, original_gff, output_faa):
    # 1. 加载“白名单”：只保留出现在过滤后 GFF 里的 Gene ID
    whitelist_gids = set()
    with open(filtered_gff, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.strip().split('\t')
            # 假设你的 mcscan.gff 第二列是 Gene_ID
            whitelist_gids.add(parts[1])

    print(f"[*] 已加载白名单基因数量: {len(whitelist_gids)}")

    # 2. 从原始 GFF 建立 Transcript_ID -> Gene_ID 映射 (仅限白名单)
    tid_to_gid = {}
    with open(original_gff, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if parts[2] in ["mRNA", "transcript"]:
                attr = parts[8]
                tid, gid = "", ""
                # 兼容多种 GFF 格式提取 ID 和 gene_id
                for a in attr.split(';'):
                    if a.startswith('ID='): tid = a.split('=')[1]
                    if a.startswith('gene_id='): gid = a.split('=')[1]

                # 关键：只有在白名单里的基因才会被映射
                if tid and gid in whitelist_gids:
                    tid_to_gid[tid] = gid

    # 3. 迭代序列，保留每个白名单 Gene_ID 最长的转录本
    gene_best_records = {} # {gid: record}

    for record in SeqIO.parse(raw_faa, "fasta"):
        tid = record.id
        gid = tid_to_gid.get(tid)

        if not gid: continue # 不在白名单里的转录本直接跳过

        # 清洗序列
        clean_seq = str(record.seq).replace('*', '').replace('.', '')
        record.seq = record.seq.__class__(clean_seq)
        seq_len = len(clean_seq)

        if gid not in gene_best_records or seq_len > len(gene_best_records[gid].seq):
            record.id = gid
            record.description = ""
            gene_best_records[gid] = record

    # 4. 输出
    with open(output_faa, "w") as out:
        SeqIO.write(gene_best_records.values(), out, "fasta")

    print(f"[*] 最终提取蛋白数量: {len(gene_best_records)}")

if __name__ == "__main__":
    # 使用: python3 filter_longest_and_rename_v3.py [原始蛋白] [过滤后的GFF] [原始GFF] [输出蛋白]
    rename_and_filter_strict(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])