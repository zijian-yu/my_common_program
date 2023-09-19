#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :run_mcscanx.py
# @Time :2023/6/3 17:14
import os
import time
import glob
import shutil
import itertools
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


class MCScanX:
    def __init__(self):
        self.gff_files = glob.glob("input/*.gff")
        self.name = self.get_name()

    # 根据gff文件，得到俩俩配对物种名字
    def get_name(self):
        sp_name = []
        for line in self.gff_files:
            sp_name.append(shutil.os.path.basename(line).split('.')[0])
        new = list(itertools.permutations(sp_name, 2))
        print("gff files cat finish!")
        return new

    # 得到 MCScanX 需要的 blast 文件
    def deal_blast(self):
        for i in range(len(self.name)):
            one, two = self.name[i]
            dfs = []
            for filename in [f"input/{one}_{one}.blast", f"input/{one}_{two}.blast", f"input/{two}_{two}.blast", f"input/{two}_{one}.blast"]:
                df = pd.read_csv(filename, sep='\t', header=None)
                dfs.append(df)
            df_concat = pd.concat(dfs, ignore_index=True)
            df_concat.to_csv(f"out_gff_blast/{one}_{two}.blast", sep='\t', header=None, index=None)
        print("Blast file cat finish!!!")

    # 将 7 列gff 文件，改为4列，并在第一列添加物种名字，再处理成 MCScanX 需要的gff文件
    def deal_gff(self):
        for line in self.gff_files:
            sp_brr = shutil.os.path.basename(line).split('.')[0]
            df = pd.read_csv(line, sep='\t', header=None)
            df = df.iloc[:, [0, 5, 1, 2]]
            df[0] = sp_brr + df[0].astype(str)
            df.to_csv(f"out_gff_blast/{sp_brr}.gff", sep='\t', header=None, index=None)
        print("gff files deal 4 col finish!")

        for one, two in self.name:
            gff_one = pd.read_csv(f"out_gff_blast/{one}.gff", sep='\t', header=None)
            gff_two = pd.read_csv(f"out_gff_blast/{two}.gff", sep='\t', header=None)
            new_gff = pd.concat([gff_one, gff_two], ignore_index=True)
            new_gff.to_csv(f"out_gff_blast/{one}_{two}.gff", sep='\t', header=None, index=None)
            print(f"cat {one}.gff {two}.gff > {one}_{two}.gff")
        print("gff files deal MCScanX finish!")

    # 跑 MCScanX
    def run_mcscanx(self):
        for one, two in self.name:
            # 使用 shutil 模块，复制、移动和删除
            print(one, two)
            shutil.copy(f"out_gff_blast/{one}_{two}.gff", ".")
            shutil.copy(f"out_gff_blast/{one}_{two}.blast", ".")

            subprocess.run(f"MCScanX {one}_{two}", shell=True)

            # synvisio-共线性要求文件（网站用）
            shutil.copy(f"{one}_{two}.collinearity", "output_synvisio")
            shutil.move(f"output_synvisio/{one}_{two}.collinearity", f"output_synvisio/{one}_{two}_collinear.collinearity")
            shutil.copy(f"{one}_{two}.gff", "output_synvisio") # rename gff file
            shutil.move(f"output_synvisio/{one}_{two}.gff", f"output_synvisio/{one}_{two}_coordinate.gff")

            shutil.move(f"{one}_{two}.tandem", "output")
            shutil.move(f"{one}_{two}.collinearity", "output")
            shutil.move(f"{one}_{two}.html", "output")

            os.remove(f"{one}_{two}.gff") # replace shutil.rmtree with os.remove
            os.remove(f"{one}_{two}.blast")


if __name__ == '__main__':
    start_time = time.time()
    mc = MCScanX()
    with ThreadPoolExecutor(max_workers=6) as executor:
        blast_future = executor.submit(mc.deal_blast)
        gff_future = executor.submit(mc.deal_gff)
        blast_future.result()
        gff_future.result()

        mc_future = executor.submit(mc.run_mcscanx)
        mc_future.result()
    end_time = time.time()
    print(f"Total running time: {end_time - start_time} seconds")

