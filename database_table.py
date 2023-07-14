#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :database_table.py
# @Time :2022/1/29 18:41
# 
# python3 database_table.py input_file output_file
import sys


file = open(sys.argv[1], 'r')
new_file = open(sys.argv[2], 'w+')
for line in file.readlines():
    new_line = line.split()
    new_file.write("<tr>")
    for hhh in new_line:

        if str(hhh) == str("."):
            new_file.write(f'''<td>{hhh}</td>''')
            # print(hhh)
        else:
            new_hhh = hhh.split('g')[0][:-2].lower()
            # new_file.write(f'''<td>{hhh}</td>''')
            # new_file.write(f'''<td><a href="static/files/{hhh}.txt"> {hhh} </a></td>''')
            new_file.write(f'''<td>{hhh}<br><a href="static/files/{new_hhh}_cds/{hhh}.txt" target="_blank">cds</a>&nbsp&nbsp&nbsp&nbsp&nbsp<a href="static/files/{new_hhh}_pep/{hhh}.txt" target="_blank">pep</a></td>''')


    new_file.write("</tr>" + '\n')








