#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :cds_num.py
# @Time :2023/6/5 16:41
# python3 cds_num.py
import os
import glob
import pandas as pd
import time


def get_id_dict(file_path):
    id_dict = {}
    with open(file_path, "r") as f:
        for line in f:
            gene_id, cds_start = line.strip().split("\t")[4:6]
            id_dict[gene_id] = cds_start
    return id_dict


def get_cds_dict(file_path, id_dict):
    cds_dict = {}
    with open(file_path, "r") as f:
        last_gene_id = "hhhh"
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "CDS":
                gene_id = fields[-1].split(";")[0].split("-")[1]  # （改） - 匹配基因id
                if gene_id in id_dict:
                    cds_dict[id_dict[gene_id]] = cds_dict.get(id_dict[gene_id], [0, 0, 0])
                    cds_dict[id_dict[gene_id]][0] += int(fields[4]) - int(fields[3]) + 1
                    cds_dict[id_dict[gene_id]][1] += 1
                    last_gene_id = gene_id
        for gene_id in cds_dict:
            cds_dict[gene_id][2] = round(cds_dict[gene_id][0]/cds_dict[gene_id][1], 2)
    return cds_dict


def read_files(file_list, cds_dict):
    df = pd.concat([pd.read_csv(file_path, sep="\t", header=None, names=[os.path.basename(file_path).split(".")[0]]) for file_path in file_list], axis=1)
    df1 = df.copy()
    df2 = df.copy()
    df3 = df.copy()
    for index, row in df.iterrows():
        for column in df.columns:
            if row[column] in cds_dict:
                df1.at[index, column] = cds_dict[row[column]][0]
                df2.at[index, column] = cds_dict[row[column]][1]
                df3.at[index, column] = cds_dict[row[column]][2]
            else:
                df1.at[index, column] = "NA"
                df2.at[index, column] = "NA"
                df3.at[index, column] = "NA"
    df1.to_csv(os.path.join("output", "cds_all_len.txt"), sep="\t", index=False)  # （可改） - 基因 cds 总长度文件
    df2.to_csv(os.path.join("output", "cds_num.txt"), sep="\t", index=False)  # （可改） - 基因 cds 数目文件
    df3.to_csv(os.path.join("output", "cds_average_len.txt"), sep="\t", index=False)  # （可改） - 基因 cds 平均长度文件


if __name__ == "__main__":
    start_time = time.time()
    id_dict = get_id_dict(os.path.join("input", "Dno.new.gff"))  # 改 - 处理之后的gff文件
    cds_dict = get_cds_dict(os.path.join("input", "old_gff.gff"), id_dict)  # 改 - 原始的gff文件
    read_files(glob.glob(os.path.join("input", "*.txt")), cds_dict)  # 改 - 四种基因的 txt 后缀文件
    end_time = time.time()
    print("The program is finished {:.2f} s.".format(end_time - start_time))
