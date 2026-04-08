import collections
import sys
# 定义一个字典来存储每个位点的转座信息
transposons = collections.defaultdict(lambda: {'te': set(), 's_count': 0, 'p_count': 0, 'reads': set()})

# 处理partial insertion 的文件
with open(sys.argv[1], 'r') as p_file:
    for line in p_file:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        te = fields[3]
        s_count = 0
        p_count = 1
        reads = "p_" + fields[4]

        # 更新位点的转座信息
        transposons[(chrom, start, end)]['te'].add(te)
        transposons[(chrom, start, end)]['p_count'] += p_count
        transposons[(chrom, start, end)]['reads'].add(reads)

# 处理spanning insertion 的文件
with open(sys.argv[2], 'r') as s_file:
    for line in s_file:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        te = fields[3]
        s_count = 1
        p_count = 0
        reads = "s_" + fields[4]

        # 更新位点的转座信息
        transposons[(chrom, start, end)]['te'].add(te)
        transposons[(chrom, start, end)]['s_count'] += s_count
        transposons[(chrom, start, end)]['reads'].add(reads)

# 将结果写入输出文件
with open(sys.argv[3], 'w') as output_file:
    for (chrom, start, end), info in transposons.items():
        te_list = ';'.join(info['te'])
        s_count = info['s_count']
        p_count = info['p_count']
        reads_list = ';'.join(info['reads'])
        output_file.write(f"{chrom}\t{start}\t{end}\t{te_list}\t{s_count}\t{p_count}\t{reads_list}\n")
