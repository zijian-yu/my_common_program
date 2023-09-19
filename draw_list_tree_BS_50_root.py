#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author : yuzijian1010@163.com
# @FileName : draw_list_tree_BS_50_root.py
# @Time : 2023/6/18 2:20
# python3 draw_list_tree_BS_50_root.py

"""
输入文件：
    tree_Pvu_Adu.txt: 共线性列表文件, 一行用于构建一个基因树
    input_pep: 对应物种的 pep 文件

输出文件：
    list_tree_BS_50_root_png.zip: 根节点图片打包
    list_tree_BS_50_root_file.zip: 根节点文件打包，用于修饰基因树

注意：
    1. 使用 muscle 和 FastTree 进行比对和构树
    2. 去掉 bootstrap 值要小于 50 的
"""
import os
import glob
import pylab
import shutil
import numpy as np
from Bio import SeqIO, Phylo
import concurrent.futures


class ListTree:
    def __init__(self, input_pep_files, vvi_file):
        self.input_pep_files = input_pep_files
        self.vvi_file = vvi_file

    def pep_files_ids(self):
        pep_gene_ids = {}
        for file in self.input_pep_files:
            file_name = os.path.basename(file).split(".")[0]
            pep_gene_ids[file_name] = {record.id: str(record.seq).replace("*", "").replace(".", "") for record in SeqIO.parse(file, "fasta")}
        return pep_gene_ids

    def write_tree_file(self, sp_geneid_dict):
        shutil.rmtree("output", ignore_errors=True)
        os.makedirs("output", exist_ok=True)
        
        with open(self.vvi_file, 'r') as mc_file:
            for line in mc_file:
                new_line = line.split()
                vvi_name = new_line[0]
                with open(f'output/{vvi_name}.txt', 'w') as new_file:
                    for hh_line in new_line:
                        if hh_line != ".":
                            gene_id = ''.join(filter(str.isalpha, hh_line.split("g")[0]))
                            new_file.write(">" + hh_line + '\n' + str(sp_geneid_dict[gene_id][hh_line]) + '\n')

    def draw_tree(self, name_file):
        """
        读取newick格式的进化树文件, 绘制进化树并保存为png格式
        :param name_file: 进化树文件名
        """
        file_name = os.path.basename(name_file).split(".")[0]
        tree = Phylo.read(f"{name_file}", "newick")
        Phylo.draw_ascii(tree)
        outgroup_clade = next((clade for clade in tree.find_clades() if clade.name == file_name), None)
        if outgroup_clade is None:
            raise ValueError(f"target '{file_name}' is not in this tree")
        tree.root_with_outgroup(outgroup_clade)  # 以外群进化枝的最近共同祖先为根节点
        Phylo.draw(tree)  # 绘制进化树
        pylab.show()  # 显示进化树
        pylab.savefig(f"{name_file}.png")  #（可改） - 保存进化树为 png 格式
        # pylab.savefig(f"{name_file}.pdf")  #（可改） - 保存进化树为 pdf 格式
        with open(f"{name_file}_rooted.newick", 'w') as rooted_file:
            Phylo.write(tree, rooted_file, 'newick')

    def run_tree_50(self, tree_name):
        os.system(f"muscle -in {tree_name} -out {tree_name}.clw")
        os.system(f"FastTree {tree_name}.clw > {tree_name}.tree")

        with open(f"{tree_name}.tree", 'r') as file_new:
            hhh_list = []
            for hhh_line in file_new.readlines():
                one = hhh_line.split(")")[1:-1]
                hhh_list.extend([float(xxx.split(":")[0]) for xxx in one if xxx.split(":")[0] != ''])
            if (np.array(hhh_list) > 0.5).all():
                with open(f'{tree_name}_50.tree', 'w') as hhh_file:
                    hhh_file.write(str(hhh_line))
                self.draw_tree(f"{tree_name}_50.tree")

    def run(self):
        sp_geneid_dict = self.pep_files_ids()
        self.write_tree_file(sp_geneid_dict)
        with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:  # (可改) - 最大进程数
            executor.map(self.run_tree_50, glob.glob(os.path.join("output", "*.txt")))


if __name__ == '__main__':
    tree = ListTree(glob.glob("input_pep/*.pep"), "tree_Pvu_Adu.txt")  # 改 - 共线性列表文件
    tree.run()

    shutil.rmtree("*.zip", ignore_errors=True)
    os.system(f"cd {os.path.join(os.getcwd(), 'output')} && zip -r ../list_tree_BS_50_root_png.zip *.png")
    os.system(f"cd {os.path.join(os.getcwd(), 'output')} && zip -r ../list_tree_BS_50_root_file.zip *rooted.newick")
