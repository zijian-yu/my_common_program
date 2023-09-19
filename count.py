#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yzj
# @FileName :run_mcscanx.py
# @Time :2023/6/3 17:34
import re
from math import *


genes = {}
genes_CDS = {}
genes_exon = {}


for row in open('Nida.genomic.gff3'):
	rows = row.strip().split()

	if rows[2] == 'gene':
		gene_id = re.split('[;=]', rows[8])[1]
		gene_length = int(rows[4]) - int(rows[3])+1
		genes[gene_id] = f'{gene_length}_{rows[6]}'

	elif rows[2] == 'CDS':
		gene_id = re.split('[;=]', rows[8])[1].split('.')[0]
		cds_length = int(rows[4]) - int(rows[3])+1
		if gene_id in genes_CDS:
			genes_CDS[gene_id].append(cds_length)
		else:
			genes_CDS[gene_id] = [cds_length]

	elif rows[2] == 'exon':
		gene_id = re.split('[;=]', rows[8])[1].split('.')[0]
		data = [int(rows[3]), int(rows[4])]
		if gene_id in genes_exon:
			genes_exon[gene_id].append(data)
		else:
			genes_exon[gene_id] = [data]


sf = open('count.csv', 'w')
sf.write(f'gene id,gene length,total exon length,exon number,max exon length,total CDS length,CDS number,max CDS length,total intron length,intron number,max intron length\n')
for n,v in genes.items():
	gene_length, chain = v.split('_')
	cds_length, cds_num, max_cds_length = sum(genes_CDS[n]), len(genes_CDS[n]), max(genes_CDS[n])

	exon_length = []
	intron_length = []
	exons = genes_exon[n]
	for ex in exons:
		exon_length.append(ex[1]-ex[0]+1)

	if chain == '-':
		for i in range(len(exons)-1):
			if len(exons) == 1: break
			exon_lens = abs(exons[i+1][1]-exons[i][0]+1)
			intron_length.append(exon_lens)

	if chain == '+':
		for i in range(len(exons)-1):
			if len(exons) == 1: break
			exon_lens = abs(exons[i+1][0]-exons[i][1]-1)
			intron_length.append(exon_lens)

	exon_length, exon_num, max_exon_length = sum(exon_length), len(exon_length), max(exon_length)
	if intron_length:
	    intron_length, intron_num, max_intron_length = sum(intron_length), len(intron_length), max(intron_length)
	else:
		intron_length, intron_num, max_intron_length = 0, 0, 0

	write_row = f'{n},{gene_length},{exon_length},{exon_num},{max_exon_length},{cds_length},{cds_num},{max_cds_length},{intron_length},{intron_num},{max_intron_length}\n'
	sf.write(write_row)


# print(genes)
# print(genes_CDS)
# print(genes_exon)