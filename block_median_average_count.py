#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@File    : block_median_count.py
@IDE     : PyCharm 
@Author  : yuzijian1010@163.com
@Date    : 2023/6/2 14:34 
@Last Date    : 2023/9/8 9:50 
"""
"""
python3 block_median_count.py
对 block、ks 文件进行处理, 得到block块的ks中位数或平均数, 然后用matlab算ks峰: ks是0-2
"""
import re
import glob
import os
import time

def colinear(ks_file, block_file, out_colinear_file, block_num, sp_list, use_median):
    with open(f"input/{ks_file}") as ks:
        hash = {line[0] + "\t" + line[1]: line[3] for line in map(str.split, ks)}

    with open(f"output/{out_colinear_file}", "w") as out, open(f"input/{block_file}") as f:
        oneblock = []
        for line_t in f:
            if line_t.startswith(("\n", "+", "overlap", "self", "the")):
                oneblock = []
                continue

            line = line_t.strip().split()
            if "".join(re.findall(r'[a-zA-Z]', str(line[0]))[:-1]) in sp_list and "".join(re.findall(r'[a-zA-Z]', str(line[2]))[:-1]) in sp_list:
                sm = line[0] + "\t" + line[2]
                if sm in hash:
                    oneblock.append(float(hash[sm]))

            elif line_t.startswith(">") and len(oneblock) >= block_num:
                if use_median:
                    array = sorted(oneblock)
                    ord = len(array) % 2
                    if ord == 0 and 2 >= float(array[len(array) // 2]) >= 0:
                        out.write(str(array[len(array) // 2]) + "\n")
                    if ord == 1:
                        od = float((float(array[(len(array) + 1) // 2]) + float(array[(len(array) - 1) // 2]))) / 2
                        if 2 >= od >= 0:
                            out.write(str(od) + "\n")
                else:
                    avg = sum(oneblock) / len(oneblock)
                    if 2 >= avg >= 0:
                        out.write(str(avg) + "\n")
                oneblock = []


if __name__ == "__main__":
    start_time = time.time()
    files = glob.glob("input/*.txt")
    name_list = list(set([os.path.basename(file).split(".")[0] for file in files]))
    sp_list = ["Nda", "Aco", "Cch", "Amtr", "Afi", "Cag", "Prh", "Epu", "Atr"]  # 改 - 物种的简称，如果不在列表中，不会进行计算
    for line in name_list:
        block_num = 8 if line.split("_")[0] == line.split("_")[1] else 6  # 改 - 目前设置种间为6，种内为8，block基因的大小
        colinear(f"{line}.ks.txt", f"{line}.block.rr.txt", f"{line}_{block_num}.txt", block_num, sp_list, use_median=True)  # 改 - use_median = True 为中位数，False 为平均数。
    end_time = time.time()
    print(f"Total running time: {end_time - start_time:.2f} seconds")

