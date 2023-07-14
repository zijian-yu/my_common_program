#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author : yuzijian1010@163.com
# @FileName : P_index.py
# @Time : 2023/6/12 15:55
# @Last Time : 2023/6/24 16:07
# python3 P_index.py
import sys
import re
import numpy as np
import glob
import os


genome_fractionation = {}
def chr_retain(arr):
    chr_ = arr[0][0]
    retain = []
    num = 0
    for i in range(0, len(arr), buchang):
        start = i - cks_2
        end = i + cks_2
        a = []
        num += 1
        for j in range(1, len(arr[i])):
            for k in range(j + 1, len(arr[i])):
                num_1 = sum([1 for index in range(start, end) if 0 <= index <= len(arr) - 1 and re.search(r"\d", arr[index][j])])
                num_2 = sum([1 for index in range(start, end) if 0 <= index <= len(arr) - 1 and re.search(r"\d", arr[index][k])])
                length = sum([1 for index in range(start, end) if 0 <= index <= len(arr) - 1])
                retain_1 = num_1 / length
                retain_2 = num_2 / length
                if retain_1 == 0 or retain_2 == 0:
                    a.append(0)
                    continue
                sigema1 = (retain_1 - retain_2) / ((retain_1 + retain_2) * 0.5)
                sigema2 = (retain_2 - retain_1) / ((retain_1 + retain_2) * 0.5)
                if sigema_start < sigema1 < sigema_end:
                    a.append(1)
                    ooo.append(sigema1)
                elif sigema_start < sigema2 < sigema_end:
                    a.append(-1)
                    uuu.append(sigema2)
                else:
                    a.append(0)
        retain.append(a)

    b = []
    genome_fractionation[arr[0][0]] = [[len(arr), num], []]
    for i in range(len(retain[0])):
        tmp = [retain[j][i] for j in range(len(retain)) if retain[j][i] != 0]  # 取出非零元素
        num1 = tmp.count(1)  # 统计1的数量
        num2 = tmp.count(-1)  # 统计-1的数量
        num0 = len(retain) - num1 - num2  # 统计0的数量
        genome_fractionation[arr[0][0]][1].extend([num1, num2, num0])  # 改 - 将结果存储到 genome_fractionation 中
        b.append([abs(sum(tmp)), len(tmp)] if tmp else [0, 0])
    return b


def calculate_retain(file_path, sub_geno_num):
    array = []
    num = 0
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if num > 0 and not line.startswith(".") and not re.search(r"\w+\d+s\d+", line):
                arr = line.split()
                ch = re.findall(r"\D+(\d*)[g|G|B]", arr[0])[0]
                array.append([ch] + arr[1:sub_geno_num+1])  # 改 - 亚基因组
            num += 1
    array = sorted(array, key=lambda x: x[0])
    all_retain = []
    list_ = []
    for i in range(len(array)):
        list_.append(array[i])
        if i == len(array) - 1 or array[i][0] != array[i + 1][0]:
            arr = chr_retain(list_)
            all_retain.append(arr)
            list_ = []
    all_ = []       
    print(f"--------------{os.path.basename(file_path).split('.')[0]}--------------")
    for i in range(len(all_retain[0])):
        num = sum([all_retain[j][i][1] for j in range(len(all_retain))])
        total = sum([all_retain[j][i][0] for j in range(len(all_retain))])
        res = total / num
        all_.append(res)
        if i < 4:  # （可能要改）
            print(f"P-index between subgenomes: {num} ----- {res}")  # 输出 - 基因数和保留度的平均值
    result = sum(all_) / len(all_)
    print(f"\nP-index: {result}\n")  # 输出 - 所有位置的保留度的平均值
    # print(f"--------------over--------------")


def Genomic_fractionation(file_path, out_file):
    with open(out_file, 'a+') as out:
        array = file_path
        array = sorted(array, key=lambda x: x[0])
        array = np.array(array)
        array = np.delete(array, 0, axis=1)
        a, b, c, d = np.sum(array, axis=0)
        arr1 = [a, b, c, d]
        for i in range(1, array.shape[1]):
            out.write(f"{arr1[i]}({arr1[i]/arr1[0]*100:.2f}%)\t")
            for j in range(array.shape[0]):
                res = array[j][i]/array[j][0]*100
                out.write(f"{array[j][i]}({res:.2f}%)\t")
            out.write('\n')
        out.write(f"{a}\t")
        for i in range(array.shape[0]):
            out.write(f"{array[i][0]}({array[i][0]/a*100:.2f}%)\t")
        out.write('\n')


if __name__ == "__main__":
    ooo = []
    uuu = []

    cks = 115  # 改 - 窗口数
    buchang = 1  # 改 - 步长值 - 依次跳几行
    # sigema_start = 0.05  # 改 - δ 的范围
    # sigema_end = 1.2  # 改 - δ 的范围
    sigema_start = 0.05  # 改 - δ 的范围
    sigema_end = 1.8  # 改 - δ 的范围
    cks_2 = cks // 2
    
    all_files = glob.glob("input_two/*")  # (可改) - 获取 input_two 文件夹中的所有文件
    for file in all_files:

        all_retain = calculate_retain(file, 2)  # 改 - 输入tab键分开的 mc 文件; 几套亚基因组！！！
        output_file = "output/" + os.path.basename(file)  # (可改) - 输出文件名字; 存储到 output 文件夹中

        os.remove(output_file) if os.path.exists(output_file) else None
        ooo = sorted(ooo, reverse=True)
        uuu = sorted(uuu, reverse=True)
        print(f"A>B Max and Min p-index:  {round(ooo[0], 5)}, {round(ooo[-1], 5)}")
        print(f"B>A Max and Min p-index:  {round(uuu[0], 5)}, {round(uuu[-1], 5)}\n")

        number = (len([j[1] for j in genome_fractionation.values()][0]))
        for j in range(0, number, 3):
            new_list = [[int(key2), value2[0][1]] + value2[1][j:j+3] for key2, value2 in genome_fractionation.items()]
            Genomic_fractionation(new_list, output_file)

