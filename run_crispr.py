#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @File    ：run_crispr.py
# @Author  ：yuzijian1010@163.com
# @Date    ：2023/6/24 17:26 
# @Last time: 2023/6/24 17:26 
# python3.6 ：run_crispr.py
# CRISPRCasFinder.pl -in Apr_genomic.fna -cas -log -out RES_Sequence -cf CasFinder-2.0.3 -def G -force -so /bio/yzj/tools/CRISPRCasFinder-release-4.2.20/sel392v2.so

# 也是可以使用 cds 来运行该流程

import os
import glob
import time
import shutil
import subprocess
import concurrent.futures


def process_species(species_name):
    cmd = f"CRISPRCasFinder.pl -in input_genome/{species_name}.fasta -cas -log -out {species_name} -cf CasFinder-2.0.3 -def G -force -so /bio/yzj/tools/CRISPRCasFinder-release-4.2.20/sel392v2.so -cpuM 12 -cpuP 12"
    subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    species_list = [os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join("input_genome", '*.fasta'))]
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:  # 改 - 最大线程数
        for species in species_list:
            shutil.copy(os.path.join("input_genome", f"{species}.fasta"), os.path.join(f"{species}.fasta"))
            executor.submit(process_species, species)
    end_time = time.time()
    print(f"The program run time is: {end_time - start_time} 秒")

