#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :draw_lollipop_chart.py
# @Time :2023/9/9 0:30
# @Last time :2023/9/9 0:30
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib


def plot_data(file_path, colors, y_range=(12, 26), output_filename='plot.png'):
    with open(file_path, 'r') as file:
        header = file.readline().strip().split('\t')
        car_name_label = header[0]
        mpg_label = header[1]
        group_label = header[2]
    
    df = pd.read_csv(file_path, sep='\t')
    grouped = df.groupby(group_label)

    fig, ax = plt.subplots()
    for i, (name, group) in enumerate(grouped):
        for x, y in zip(group[car_name_label], group[mpg_label]):
            ax.plot([x, x], [y, 0], color='black', linestyle='-', linewidth=1.5, zorder=1)
            ax.annotate(f'{y}', (x, y), textcoords="offset points", xytext=(0, -2), ha='center', zorder=3, fontsize=8)
        ax.scatter(group[car_name_label], group[mpg_label], marker='o', s=250, label=f'Group {name}', color=colors[i], zorder=2)

    # ax.set_ylim(12, 26)  # 改 - 设置y轴的起始和终止坐标
    ax.set_ylim(*y_range)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(y=16, color='gray', linestyle='--', linewidth=1)  # 改 - 灰色线的位置
    ax.axhline(y=22, color='gray', linestyle='--', linewidth=1)  # 改 - 灰色线的位置
    plt.xticks(rotation=45, ha="right")
    plt.subplots_adjust(right=0.8)  # 调整子图布局，将标签和标题移到图片框的右侧以外
    # 使用从文件中获取的标题名字
    ax.set_xlabel(car_name_label)
    ax.set_ylabel(mpg_label)
    ax.set_title(f'{mpg_label} by {car_name_label} and {group_label}')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig('plot.png', dpi=300)


if __name__ == "__main__":
    matplotlib.rcParams['font.family'] = 'Arial'  # 改 - 设置字体为 Times New Roman,Arial
    colors = ['steelblue', 'limegreen', 'tomato', 'gold']  # 改 - 自定义颜色列表，每个颜色对应一个组
    plot_data('test.txt', colors, y_range=(12, 26), output_filename='plot.png')  # 改 - 1.文件名字. 2.定义的颜色. 3.y轴范围. 4. 输出文件名字
