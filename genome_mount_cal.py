#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :genome_mount_cal.py
# @Time :2023/8/30 23:15
# python3 genome_mount_cal.py
# 计算基因组的挂载率，删除了小于2000bp的contig
from Bio import SeqIO


def calculate_gene_length_ratio(fasta_file, num_chromosomes):
    sequences = SeqIO.parse(fasta_file, "fasta")
    gene_lengths = {}
    for record in sequences:
        gene_id = record.id
        gene_length = len(record.seq)
        if gene_length > 2000:  # （改）- 注意，删除了小于2000bp的contig
            gene_lengths[gene_id] = gene_length
    # 按碱基数目从大到小对基因进行排序
    sorted_genes = sorted(gene_lengths.items(), key=lambda x: x[1], reverse=True)
    # 取前12个基因的碱基数目并计算总长度
    top_genes = sorted_genes[:num_chromosomes]
    total_top_length = sum([length for _, length in top_genes])
    total_genome_length = sum(gene_lengths.values())
    ratio = total_top_length / total_genome_length
    return ratio


if __name__ == "__main__":
    fasta_file = "contig_HiC.fasta"  # 改 - 基因组文件名字
    num_chromosomes = 12  # 染色体数目
    result = calculate_gene_length_ratio(fasta_file, num_chromosomes)
    print(f"Top {num_chromosomes} gene length ratio: {result:.4f}")

