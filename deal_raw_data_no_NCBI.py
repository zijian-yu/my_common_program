#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@File    ：deal_raw_data_no_NCBI.py
@IDE     ：PyCharm 
@Author  ：YZJ
@Date    ：2023/6/1 6:13 
"""
"""
python3 deal_raw_data_no_NCBI.py
"""
import sys
import re
from Bio import SeqIO


def rewrite_gff(gff_file, output_file, feature_type, name_tag):
    with open(gff_file, 'r') as IN, open(output_file + '.new.gff', 'w') as OUT, open(output_file + ".lens", 'w') as lens_file:
        hash = {}
        num = 0
        old = ''
        gene_to_gene = {}
        total_genes = 0

        gff_dish = {}
        lens_dish = {}
        for line in IN:
            if line.startswith('#'):
                continue
            arr = line.split()
            if feature_type in arr[2]:
                a = arr[8].split(';')
                a[1] = a[1].replace(name_tag, '')
                if a[1] in hash:
                    continue
                hash[a[1]] = 1
                arr0_num = ''.join(re.findall(r'\d+', arr[0]))
                chr_input = output_file + arr0_num
                num = num + 1 if old == chr_input else 1
                lens_dish[chr_input] = num
                gff_dish[f"{chr_input}g{num:04d}"] = f"{chr_input}\t{arr[3]}\t{arr[4]}\t{arr[6]}\t{a[1]}\t{chr_input}g{num:04d}\t{num}\n"
                old = chr_input
                gene_to_gene[a[1]] = f"{chr_input}g{num:04d}"
                total_genes += 1

        gff_dish = sorted(gff_dish.items(), key=lambda x: x[0]) # 根据基因id对gff排序
        for _, line in gff_dish:
            OUT.write(line)

        lens_dish = sorted(lens_dish.items()) # 根据染色体号排序
        qian_num = 200
        for line, lens in lens_dish:
            lens_file.write(f"{line}\t{lens}\t{round(lens / total_genes * 2000 + qian_num)}\n")
            qian_num = round(lens / total_genes * 2000 + qian_num)  # 四舍五入
        print("gff and lens files is run over!!!")
    return gene_to_gene


def replace_gene_id(cds_file, pep_file, sp_name, gene_to_gene):
    cds_records = SeqIO.parse(cds_file, "fasta")
    pep_records = SeqIO.parse(pep_file, "fasta")
    cds_new_records = []
    pep_new_records = []
    for cds_record, pep_record in zip(cds_records, pep_records):
        gene_id = cds_record.id.split('.')[0]
        if gene_id in gene_to_gene:
            cds_record.id = gene_to_gene[gene_id]
            cds_record.description = ''
            pep_record.id = gene_to_gene[gene_id]
            pep_record.description = ''
            cds_new_records.append(cds_record)
            pep_new_records.append(pep_record)
    SeqIO.write(cds_new_records, sp_name + '.cds', "fasta")
    SeqIO.write(pep_new_records, sp_name + '.pep', "fasta")
    print("cds and pep files is replace over!!!")


if __name__ == "__main__":
    hh = rewrite_gff("old_gff.gff3", "Lpu", 'gene', 'Name=')  # 改 - 原始gff文件，物种简称，gene or mRNA，匹配gff原始基因id
    replace_gene_id("old_cds.fasta", "old_pep.fasta","Lpu", hh)  # 改 - 原始cds文件，原始pep文件，物种简称，对应基因id关系

