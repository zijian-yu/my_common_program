#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_interproscan.py
# @Time :2023/6/3 0:29
# python3 run_interproscan.py
import glob
import subprocess
import concurrent.futures
import time


def run_interproscan(input_file, CPU):
    cmd = ["interproscan.sh", "-i", input_file, "-b", "result/", "-goterms", "-iprlookup", "-pa", "-dp", "-f", "tsv", "-cpu", str(CPU)]
    print("Running command: ", " ".join(cmd))
    subprocess.run(cmd)


if __name__ == '__main__':
    start_time = time.time()
    pep_file = glob.glob("input/*.pep")  # （可改） - 文件的路径和文件后缀
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:  # 改 - 线程池最大线程数
        for file in pep_file:
            executor.submit(run_interproscan, file, 10)  # 改 - 运行单个文件所占的CPU
    end_time = time.time()
    print(f"Time taken: {end_time - start_time} seconds")
