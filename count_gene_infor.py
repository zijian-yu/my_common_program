#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@File    : count_gene_infor.py
@IDE     : PyCharm 
@Author  : YZJ
@Date    : 2023/5/31 10:43 
"""
"""
python3 count_gene_infor.py
"""
import glob
import os


def get_counts(fn_gff):
    counts = {"gn": 0, "gs": 0, "en": 0, "es": 0, "un": 0, "us": 0, "mn": 0, "ms": 0, "in_": 0, "is_": 0, "cn": 0, "cs": 0}
    deal = ("#", " ", "\n")
    with open(fn_gff) as f:
        for line in f:
            if line.startswith(deal):
                continue
            a = line.split()
            # print(a[2])
            if a[2] == "gene":
                counts["gn"] += 1
                counts["gs"] += abs(int(a[4]) - int(a[3])) + 1
            elif a[2] == "mRNA":
                counts["mn"] += 1
                counts["ms"] += abs(int(a[4]) - int(a[3])) + 1
            elif a[2] == "exon":
                counts["en"] += 1
                counts["es"] += abs(int(a[4]) - int(a[3])) + 1
            elif a[2] == "CDS":
                counts["cn"] += 1
                counts["cs"] += abs(int(a[4]) - int(a[3])) + 1
            elif "UTR" in a[2]:
                counts["un"] += 1
                counts["us"] += abs(int(a[4]) - int(a[3])) + 1
    if counts["en"] == 0:
        with open(fn_gff) as f:
            ls = 0
            le = 0
            for line in f:
                if line.startswith(deal):
                    continue
                a = line.split()
                if "UTR" in a[2] or a[2] == "CDS":
                    counts["en"] += 1
                    counts["es"] += abs(int(a[4]) - int(a[3])) + 1
                    if abs(ls - int(a[4])) == 1 or abs(le - int(a[3])) == 1:
                        counts["en"] -= 1
                    ls = int(a[3])
                    le = int(a[4])
    counts["in_"] = counts["en"] - counts["mn"]
    counts["is_"] = counts["ms"] - counts["es"]
    return tuple(counts.values())


if __name__ == "__main__":
    files = glob.glob("input/*.gff")  # （可改） - 输入文件夹的名字

    with open("counts.txt", "w") as f:  # （可改） - 最终统计文件名字
        f.write("filename\tgene_num\tgene_len\texon_num\texon_len\tUTR_num\tUTR_len\tmRNA_num\tmRNA_len\texon_num\texon_len\tCDS_num\tCDS_len\n")
        for fn_gff in files:
            print(fn_gff) 
            filename = os.path.basename(fn_gff)
            counts = get_counts(fn_gff)
            f.write(f"{filename}\t{counts[0]}\t{counts[1]}\t{counts[2]}\t{counts[3]}\t{counts[4]}\t{counts[5]}\t{counts[6]}\t{counts[7]}\t{counts[8]}\t{counts[9]}\t{counts[10]}\t{counts[11]}\n")
