#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :GeneCollinearity_gene.py
# @Time :2023/6/15 11:21
# python3 GeneCollinearity_gene.py

"""
1. 基因id文件中, ">"号可有可无。
2. 处理的是 MCScanX 的结果文件。
3. 功能: 统计基因id文件中的基因id的基因对数目; 所有的基因id数目; 前两的比值; 基因对的平均e值
"""
import pandas as pd


class GeneCollinearity:
    def __init__(self, file_path, gene_file_path, double_single=True):
        self.file_path = file_path
        self.gene_list = self.get_gene_list(gene_file_path)
        self.double_single = double_single
        self.all_geneid = len(self.gene_list)

    def get_gene_list(self, file_path):
        with open(file_path, 'r') as f:
            gene_list = [line.strip() for line in f.readlines()]
            gene_list = pd.read_csv(file_path, header=None, sep='>')[1].str.strip().tolist() if ">" in gene_list[0] else [gene.strip() for gene in gene_list]
        return gene_list

    def get_collinearity_dict(self):
        e_value_list = []
        e_num = 0
        with open(self.file_path, 'r') as f:
            for line in f.readlines():
                if line.startswith('##'):
                    continue
                line = line.strip().split('\t')
                gene1, gene2 = line[1:3]
                if (self.double_single and gene1 in self.gene_list and gene2 in self.gene_list) or (not self.double_single and (gene1 in self.gene_list or gene2 in self.gene_list)):
                    e_num += 1
                    e_value_list.append(float(line[3]))
        avg_e_value = sum(e_value_list) / len(e_value_list)
        all_geneid_num = self.all_geneid
        return e_num, avg_e_value, all_geneid_num


if __name__ == "__main__":
    # 改 - MCScanX结果文件; 基因id文件; double_single=True:两个基因都在，double_single=False:一个基因在。
    gene_collinearity = GeneCollinearity('Pv_Pv.collinearity', 'Pv.gene.txt', double_single=True)
    e_num, avg_e_value, all_geneid_num = gene_collinearity.get_collinearity_dict()

    print("Gene id is the number of genes in the block:", e_num)  # 基因id在block中的数目
    print("Count all the genes in the block file:", all_geneid_num)  # 基因列表中的所有基因id数目
    print("Gene id is the percentage of the number of genes in the block (%):", round(e_num / all_geneid_num * 100, 5))  # 前两比值
    print("The average e-value of the gene id in the block:", avg_e_value)  # 平均的e值

