#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :piliang_quchong.py
# @Time :2022/4/11 0:40
import os
from Bio import SeqIO

path = os.getcwd()
files = os.listdir(f"""{path}/files""")
os.mkdir("new_files")
for line in files:
    new_line = line.split(".")[0]
    new_file = open(f"""new_files/{new_line}.new.pep""", 'w')
    for sequence in SeqIO.parse(f"""files/{new_line}.pep""", 'fasta'):
        if str(sequence.id) == str("000"):
            continue
        else:
            new_file.write('>' + str(sequence.id) + '\n' + str(sequence.seq) + '\n')
