#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :GeneCollinearity_block.py
# @Time :2023/6/15 11:21
# python3 GeneCollinearity_block.py

"""
1. 基因id文件中, ">"号可有可无。
2. 处理的是 MCScanX 的结果文件。
3. 功能: 统计基因id文件中的基因id所在block块中的所有基因; block块中的所有基因; 前两的比值; 所在block块的平均e值
"""
import pandas as pd


class GeneCollinearity:
    def __init__(self, file_path, gene_file_path):
        self.file_path = file_path
        self.gene_file_path = gene_file_path

        with open(self.gene_file_path, 'r') as f:  # 处理基因id文件中有无">"号不
            gene_list = [line.strip() for line in f.readlines()]
            gene_list = pd.read_csv(self.gene_file_path, header=None, sep='>')[1].str.strip().tolist() if ">" in gene_list[0] else [gene.strip() for gene in gene_list]
        self.gene_list = gene_list

        self.collinearity_dict, self.e_value_dict, self.jishu_dict, self.all_gene_id = self.get_collinearity_dict()

    def get_collinearity_dict(self):
        collinearity_dict, e_value_dict, e_num, jishu_dict, all_gene_id = {}, {}, 1, {}, []
        with open(self.file_path, 'r') as f:
            for line in f.readlines():
                if line.startswith('##'):
                    e_value_dict[e_num], e_num, collinearity_dict[e_num] = line.split()[4].split("=")[1], e_num+1, []
                    continue
                line = line.strip().split('\t')
                collinearity_dict[e_num].extend(line[1:3])
                all_gene_id.extend(line[1:3])
                if all(gene in self.gene_list for gene in line[1:3]):
                    jishu_dict.setdefault(e_num, []).extend(line[1:3])
        return collinearity_dict, e_value_dict, jishu_dict, all_gene_id

    def get_pp2a_coll_list(self):
        return [gene_id for key in self.jishu_dict for gene_id in self.collinearity_dict[key]]

    def get_gene_id_count(self):
        return len(set(self.get_pp2a_coll_list()))

    def get_all_gene_id(self):
        return len(set(self.all_gene_id))

    def get_gene_id_percentage(self):
        return round(self.get_gene_id_count() / self.get_all_gene_id() * 100, 5)

    def get_e_values(self):
        return [float(self.e_value_dict[key]) for key in self.jishu_dict]

    def get_avg_e_value(self):
        return sum(self.get_e_values()) / len(self.get_e_values())


if __name__ == "__main__":
    gene_collinearity = GeneCollinearity('Pv_Pv.collinearity', 'Pv.gene.txt')  # 改 - MCScanX结果文件，基因id文件
    gene_id_count = gene_collinearity.get_gene_id_count()
    all_gene_id = gene_collinearity.get_all_gene_id()
    gene_id_percentage = gene_collinearity.get_gene_id_percentage()
    avg_e_value = gene_collinearity.get_avg_e_value()

    print("Gene id is the number of genes in the block:", gene_id_count)
    print("Count all the genes in the block file:", all_gene_id)
    print("Gene id is the percentage of the number of genes in the block (%):", gene_id_percentage)
    print("The average e-value of the gene id in the block:", avg_e_value)

