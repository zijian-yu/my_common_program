#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :draw_box.py
# @Time :2023/6/5 16:41
# python3 draw_box.py
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import glob
import time


def draw_box(input_file, output_file):
    sns.set(style="whitegrid")
    df1 = pd.read_csv(input_file, sep="\t")
    df1 = df1.dropna()
    df1.columns = ["Core", "Soft core", "Dispensable", "Private"]
    df1 = np.log2(df1)
    colors = ["#F9766C", "#7CAF00", "#00BFC5", "#C77CFF"]  # 改 - 颜色
    plt.figure(figsize=(8, 6), dpi=300)
    sns.boxplot(data=df1.melt(var_name='columns', value_name='values'), x='columns', y='values', 
                flierprops=dict(markerfacecolor='none', markeredgecolor='pink', marker='o', markersize=5), palette=colors, 
                width=0.5)
    plt.xticks(rotation=0)
    plt.ylabel("CDS Lengths (x1000)")
    plt.ylim(1, 100)
    plt.yticks(np.arange(1, 100, 1))
    plt.gca().autoscale(enable=True, axis='y')

    plt.gca().xaxis.set_tick_params(length=5, width=1, which='both', direction='out', bottom=True, top=False, labelbottom=True, color='darkgray')
    plt.gca().yaxis.set_tick_params(length=5, width=1, which='both', direction='out', left=True, right=False, labelleft=True, color='darkgray')
    plt.gca().yaxis.grid(True, linestyle='--', linewidth=0.5, color='darkgray', alpha=0)
    plt.savefig(output_file + ".png")
    plt.savefig(output_file + ".pdf")


if __name__ == "__main__":
    start_time = time.time()
    input_file = glob.glob(os.path.join("output", "*.txt"))
    input_file = [f for f in input_file if f.endswith(".txt")]
    for file in input_file:
        output_file = os.path.splitext(file)[0]
        draw_box(file, output_file)
    end_time = time.time()
    print("The program is finished {:.2f} s.".format(end_time - start_time))
