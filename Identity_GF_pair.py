#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :identity_GF_pair.py
# @Time :2023/6/12 16:41
# @Last time: 2023/6/23 16:38
# python3 Identity_GF_pair.py
import re
import os
import glob
import time
import shutil
import subprocess
import pandas as pd
import concurrent.futures
import matplotlib.pyplot as plt
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# 1 批量执行 blast
class BlastRunner:
    def __init__(self, max_workers):
        self.max_workers = max_workers

    def copy_input_files(self):
        input_files = glob.glob(os.path.join("input_pep_file", "*"))
        for file in input_files:
            if os.path.exists(os.path.basename(file)):
                continue
            shutil.copy(file, os.getcwd())

    def create_database(self, database_name, fasta_file):
        subprocess.run(["makeblastdb", "-in", fasta_file, "-dbtype", "prot", "-out", database_name])

    def run_blast(self, query_file, database_name, output_file):
        subprocess.run(["blastp", "-query", query_file, "-db", database_name, "-out", output_file, "-outfmt", "6", "-evalue", "1e-3", "-num_threads", "8"])  # ()

    def delete_files(self):
        # 删除除了 .py 和 .blast 文件之外的所有文件，将 .blast 文件移动到 out_blast 文件夹中
        files = glob.glob("*")
        for file in files:
            if not file.endswith(".py") and not file.endswith(".blast"):
                if os.path.isfile(file):
                    os.remove(file)
        os.makedirs("output_1_blast", exist_ok=True)
        blast_files = glob.glob("*.blast")
        for file in blast_files:
            shutil.move(file, os.path.join("output_1_blast", file))

    def main_run_blast(self, query_file, fasta_files):
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = []
            for fasta_file in fasta_files:
                futures.append(executor.submit(self.create_database, fasta_file.split(".")[0], fasta_file))
            concurrent.futures.wait(futures)

            futures = []
            for query in query_file:
                for fasta_file in fasta_files:
                    futures.append(executor.submit(self.run_blast, query, fasta_file.split(".")[0], f"{query.split('.')[0]}_{fasta_file.split('.')[0]}.blast"))
            concurrent.futures.wait(futures)
        self.delete_files()


# 2 和 3 提取基因id蛋白文件。结果写入到 output_2_pep_id 文件夹中
class GeneIdExtractor:
    def __init__(self, input_pep_files, blast_files):
        self.input_pep_files = input_pep_files
        self.blast_files = blast_files
        self.pep_id = self.pep_files_ids()

    def pep_files_ids(self):
        gene_ids = {}
        for file in self.input_pep_files:
            file_name = os.path.basename(file).split(".")[0]
            gene_ids[file_name] = {record.id: str(record.seq).replace("*", "").replace(".", "") for record in SeqIO.parse(file, "fasta")}
        return gene_ids

    def extract_gene_ids(self):
        os.makedirs("output_2_pep_id", exist_ok=True)
        for file in self.blast_files:
            file_name = os.path.basename(file).split(".")[0]
            df = pd.read_csv(file, sep='\t', header=None)
            gene_list = list(set(df[1]))  # pep_files_ids
            gene_id_file = os.path.join("output_2_pep_id", file_name + ".fasta")
            if os.path.exists(gene_id_file):
                os.remove(gene_id_file)
            for gene in gene_list:
                if file_name.split("_")[1] in self.pep_id.keys() and gene in self.pep_id[file_name.split("_")[1]].keys():
                    with open(gene_id_file, "a") as f:
                        record = SeqRecord(Seq(self.pep_id[file_name.split("_")[1]][gene]), id=gene, description="")
                        SeqIO.write(record, f, "fasta")


# 4 interproscan
class InterproscanRunner:
    def __init__(self, input_file, CPU):
        self.input_file = input_file
        self.CPU = CPU

    def run_interproscan(self):
        cmd = ["interproscan.sh", "-i", self.input_file, "-b", "output_3_interproscan/", "-goterms", "-iprlookup", "-pa", "-dp", "-f", "tsv", "-cpu", str(self.CPU)]
        subprocess.run(cmd)


# 5、6、7：筛选、提取 pep文件、统计个数
class Deal_Inter:
    def __init__(self, input_pep_files, input_cds_files, file_path, pfam_num, interpro_type):
        self.input_pep_files = input_pep_files
        self.input_cds_files = input_cds_files
        self.file_path = file_path
        self.pfam_num = pfam_num
        self.pep_gene_ids = self.pep_files_ids()
        self.cds_gene_ids = self.cds_files_ids()
        self.interpro_type = interpro_type

    def pep_files_ids(self):
        pep_gene_ids = {}
        for file in self.input_pep_files:
            file_name = os.path.basename(file).split(".")[0]
            pep_gene_ids[file_name] = {record.id: str(record.seq).replace("*", "").replace(".", "") for record in SeqIO.parse(file, "fasta")}
        return pep_gene_ids
    
    def cds_files_ids(self):
        cds_gene_ids = {}
        for file in self.input_cds_files:
            file_name = os.path.basename(file).split(".")[0]
            cds_gene_ids[file_name] = {record.id: str(record.seq).replace("*", "").replace(".", "") for record in SeqIO.parse(file, "fasta")}
        return cds_gene_ids
    
    
    def process_file(self):
        file_name = os.path.basename(self.file_path)

        with open(self.file_path, 'r') as inter_file, open(f"output_4_deal_inter/{file_name}", 'w') as new_file, open(f"output_5_inter_pep/{file_name.split('.')[0]}.fasta", 'w') as new_pep, open(f"output_5_inter_cds/{file_name.split('.')[0]}.fasta", 'w') as new_cds, open("interproscan_pep_count.txt", 'a') as interproscan_pep_count:
            number = 0
            gene_id = []
            for inter_line in inter_file:
                new_line = inter_line.strip().split("\t")
                if new_line[3] == self.interpro_type and new_line[4] in self.pfam_num[file_name.split(".")[0].split("_")[0]]:
                    new_file.write("\t".join(new_line) + "\n")
                    gene_id_list = self.pep_gene_ids[file_name.split(".")[0].split("_")[1]].keys()
                    if new_line[0] in gene_id_list and new_line[0] not in gene_id:
                        new_pep.write(f">{new_line[0]}\n{self.pep_gene_ids[file_name.split('.')[0].split('_')[1]][new_line[0]]}\n")
                        new_cds.write(f">{new_line[0]}\n{self.cds_gene_ids[file_name.split('.')[0].split('_')[1]][new_line[0]]}\n")
                        number += 1
                        gene_id.append(new_line[0])
            interproscan_pep_count.write(f"{file_name.split('.')[0]}\t{number}\n")


# 8 构建基因树文件，并画图
class GeneTree:
    def __init__(self, output_inter_5):
        self.output_inter_5 = output_inter_5

    def merge_files(self, file_name):
        # 合并所有以 “.fasta” 为后缀的文件
        fasta_files = glob.glob(f"{file_name}*.fasta")
        with open(file_name + ".fasta", "w") as outfile:
            for file in fasta_files:
                outfile.write(open(file).read())
                if file != file_name:
                    os.remove(file)

    def copy_input_files(self):
        # 复制 input_pep_file 中的文件，到本路径
        file_name = []
        for file in self.output_inter_5:
            if os.path.exists(os.path.basename(file)):
                continue
            shutil.copy(file, os.getcwd())
            file_name.append(os.path.basename(file).split("_")[0])
        for line in list(set(file_name)):
            print(line)
            self.merge_files(line)

    def draw_tree(self, tree_file):
        tree = Phylo.read(tree_file, 'newick')
        plt.rcParams['figure.figsize'] = (20, 10)
        Phylo.draw(tree, label_func=lambda x: x.name, show_confidence=False, axes=plt.gca())
        plt.savefig(tree_file.split(".")[0] + ".png")
        plt.close()

    def creat_tree(self, protein_file):
        # .fas: 序列比对文件；.1.fas：去gap文件；.contree 和 iqtree：基因树文件
        shutil.rmtree("output_6_seq_treeimg", ignore_errors=True)
        os.makedirs("output_6_seq_treeimg")
        os.system("mafft --auto " + protein_file + " > " + protein_file.split(".")[0] + ".fas")
        os.system("trimal -in " + protein_file.split(".")[0] + ".fas" + " -out " + protein_file.split(".")[0] + ".1.fas" + " -gt " + "0.3")
        os.system("iqtree2 -s " + protein_file.split(".")[0] + ".1.fas" + " -m " + "MFP" + " -bb " + "1000" + " -T" + " 8")
        self.draw_tree(protein_file.split(".")[0] + ".1.fas.contree")

    def move_files(self):
        for file in os.listdir(os.getcwd()):
            if file.endswith((".fas", ".1.fas", ".contree", ".iqtree", ".fasta", ".png")):
                shutil.move(file, "output_6_seq_treeimg")
            elif os.path.isfile(file) and not file.endswith((".py", ".txt")):
                os.remove(file)

    def run(self):
        self.copy_input_files()
        with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:  # 改 - 设置最大线程数
            for file in glob.glob("*.fasta"):
                print(file)
                gene_ids = [record.id for record in SeqIO.parse(file, "fasta")]  # 改 - 使用列表推导式
                if len(gene_ids) <= 200:  # 改 - 判断基因id数目是否小于等于 200
                    self.creat_tree(file)
        self.move_files()


# 9 读取 block 文件，提取家族相关基因对
class GeneProcessor:
    def __init__(self, fasta_files, block_files):
        self.fasta_files = fasta_files
        self.block_files = block_files
        self.genes_dict = self.get_gene_ids()
        self.move_files()

    def move_files(self):
        for file in self.fasta_files:
            shutil.copy(file, os.getcwd())

    def get_gene_ids(self):
        gene_ids = {}
        for file in self.fasta_files:
            gene_ids[os.path.basename(file.split(".")[0])] = []
            for record in SeqIO.parse(file, "fasta"):
                gene_ids[os.path.basename(file.split(".")[0])].append(record.id)
        return gene_ids

    def process_colinearscan_results(self, file_name):
        shutil.copy(os.path.join("input_block_file", os.path.basename(file_name)), os.path.join(os.getcwd(), os.path.basename(file_name)), follow_symlinks=True)

        file_name = os.path.basename(file_name)
        with open(file_name, "r") as file:
            deal = ("over", "+++", "\n", "the", " ", "self")
            block_dict = {}
            for line in file:
                if line.startswith(deal):
                    continue
                if line.startswith(">"):
                    for gene_fam, genes_id in self.genes_dict.items():
                        gene_pei_dui = []
                        gene_fam_name = gene_fam
                        for genes_id in genes_id:
                            sp_chr, sp_chr_number = genes_id.split("g")
                            if sp_chr in block_dict and int(block_dict[sp_chr][0]) > int(sp_chr_number) > int(block_dict[sp_chr][1]):
                                gene_pei_dui.append(genes_id)
                        for i in range(len(gene_pei_dui)):
                            for j in range(i+1, len(gene_pei_dui)):
                                with open(os.path.join("output_7_pair", f"{gene_fam}_{file_name.split('.')[0]}.txt"), "a+") as f:
                                    f.write(f"{gene_pei_dui[i]}\t{gene_pei_dui[j]}\n")
                    block_dict = {}
                    continue

                new_line = line.strip().split()
                gene1, gene2 = new_line[0], new_line[2]
                for key, value in self.genes_dict.items(): # 遍历genes_dict中的key和value
                    if gene1 in value and gene2 in value: # 如果gene1和gene2都在value中
                        with open(os.path.join("output_7_pair", f"{key}_{file_name.split('.')[0]}.txt"), "a+") as f: # 打开文件
                            f.write(f"{gene1}\t{gene2}\n") # 写入gene1和gene2
                for gene in [gene1, gene2]: # 遍历gene1和gene2
                    gene_key, gene_value = gene.split("g")[0], [int(gene.split("g")[1]), int(gene.split("g")[1])] # 将gene分割成key和value
                    block_dict[gene_key] = [max(gene_value[0], block_dict.get(gene_key, [0, 0])[0]), min(gene_value[1], block_dict.get(gene_key, [float('inf'), float('inf')])[1])] # 将gene_key和gene_value添加到block_dict中，并更新gene_value[0]和gene_value[1]
        os.remove(file_name)

    def run(self):
        shutil.rmtree("output_7_pair", ignore_errors=True)
        os.makedirs("output_7_pair", exist_ok=True)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for file in self.block_files:
                executor.submit(self.process_colinearscan_results, file)


# 10. 对去基因 pairs 进行去重，并合并文件  
class GFBlockAll:
    def __init__(self, pair_files):
        self.pair_files = pair_files
        self.output_dir = "output_7_pair_result"
        self.gf_name = None

    def process_file(self, line):
        data = {}
        for hhh in line:
            with open(hhh, "r") as file_open:
                self.gf_name = os.path.basename(hhh).split("_")[0]
                data[self.gf_name] = []

        for hhh in line:
            with open(hhh, "r") as file_open:
                self.gf_name = os.path.basename(hhh).split("_")[0]
                for line in file_open:
                    a, b = line.split()
                    if a != b:
                        data[self.gf_name].append(f"{a}_{b}")

        for key, value in data.items():
            with open(os.path.join(self.output_dir, f"{key}_pair.txt"), "a") as cct_gene:
                for pair in set(value):
                    a, b = pair.split("_")
                    if a != b and f"{b}_{a}" not in set(value):
                        cct_gene.write(f"{a}\t{b}\n")

    def run(self):
        shutil.rmtree(self.output_dir, ignore_errors=True)
        os.makedirs(self.output_dir, exist_ok=True)
        self.process_file(self.pair_files)


# 11 读取 input_gff_file 文件夹中的gff文件，并写入一个字典中，要求文件名字为键，嵌套一个字典，嵌套字典的键为gff的倒数第二列，嵌套字典的值为第二列和第三列，分别为gff文件的起始和终止位置
class GeneFinder:
    def __init__(self, gff_files, seq_treeimg_files):
        self.gff_files = gff_files
        self.seq_treeimg_files = seq_treeimg_files
        self.output_dir = "output_8_tandem"

    def move_files(self, files):
        for file in files:
            shutil.copy(file, os.getcwd())

    def read_gff_files(self):
        gff_dict = {}
        for file in self.gff_files:
            filename = os.path.basename(file).split(".")[0]
            with open(file, "r") as f:
                gff_dict[filename] = {}
                for line in map(str.strip, f):
                    line = line.split("\t")
                    gff_dict[filename][line[-2]] = (line[1], line[2])
        return gff_dict

    def get_gene_ids(self, seq_treeimg, gff_dict):
        gf_name = os.path.basename(seq_treeimg.split(".")[0])
        with open(seq_treeimg, "r") as f:
            for line in f:
                gene1, gene2 = line.strip().split("\t")
                sp_chr1, sp_chr2 = gene1.split("g")[0], gene2.split("g")[0]
                if sp_chr1 == sp_chr2:
                    if gene1 > gene2:
                        gene1, gene2 = gene2, gene1
                    sp_name = re.findall(r'[a-zA-Z]+', sp_chr1)[0]
                    xiao = gff_dict[sp_name][gene1][1]
                    da = gff_dict[sp_name][gene2][0]
                    chazhi = int(da) - int(xiao)
                    if chazhi <= 50000:  # 改 - 50kb以内的基因定义为串联重复基因
                        with open(f"output_8_tandem/{gf_name}.txt", "a") as f:
                            f.write(f"{gene1}\t{gene2}\n")

    def remove_block_rr_files(self):
        for file in glob.glob("*.block.rr.txt"):
            if file.endswith(".block.rr.txt"):
                os.remove(file)

    def run(self):
        self.move_files(glob.glob("output_7_pair_result/*.txt"))  # 移动文件
        gff_dict = self.read_gff_files()  # {物种：{id：seq}}

        shutil.rmtree(self.output_dir, ignore_errors=True)
        os.makedirs(self.output_dir, exist_ok=True)

        with concurrent.futures.ProcessPoolExecutor() as executor:
            for i in glob.glob("*_pair.txt"):
                executor.submit(self.get_gene_ids, i, gff_dict)  # 基因对列表

        self.remove_block_rr_files()


if __name__ == "__main__":
    """
    input_block_file: 输入 block 文件
    input_gff_file: 输入 gff 文件
    input_pep_file: 输入 物种蛋白文件(.pep)和基因家族文件(.fasta)

    output_1_blast: 基因家族和物种蛋白, blastp 比对文件
    output_2_pep_id: 基因家族 blastp 提取的 pep 序列
    output_3_interproscan: interproscan 运行结果
    output_4_deal_inter: 筛选 interproscan 结果
    output_5_inter_pep: GF_SP.fasta, interproscan 蛋白序列
    output_6_seq_treeimg: Iqtree等构树结果, 以及基因树图
    output_7_pair: 从 block 文件提取的基因对结果文件
    output_7_pair_result: 筛选去重后的结果
    output_8_tandem: tandem 基因对
    """
    start_time_all = time.time()

    # 1 批量执行 blast
    blast_runner = BlastRunner(max_workers=6)  # （可改） - 线程池最大线程数
    blast_runner.copy_input_files()
    blast_runner.main_run_blast(glob.glob("*.fasta") , glob.glob("*.pep"))  # （可改） - 家族蛋白文件, 物种蛋白文件

    # 2 和 3 提取基因id蛋白文件
    run_two_and_three = GeneIdExtractor(glob.glob("input_pep_file/*"), glob.glob("output_1_blast/*.blast"))
    run_two_and_three.extract_gene_ids()

    # 4 interproscan
    start_time = time.time()
    shutil.rmtree("output_3_interproscan", ignore_errors=True)
    os.makedirs("output_3_interproscan")
    pep_file = glob.glob("output_2_pep_id/*.fasta")  # （可改） - 文件的路径和文件后缀
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:  # 改 - 线程池最大线程数
        for file in pep_file:
            interproscan_runner = InterproscanRunner(file, 6)  # 注意：最大线程数与这个数字相乘不建议超过服务器线程数！
            executor.submit(interproscan_runner.run_interproscan)
    end_time = time.time()
    print(f"Interproscan time taken: {end_time - start_time} seconds")

    # 5、6、7：筛选、提取 pep文件、统计个数
    shutil.rmtree("output_4_deal_inter", ignore_errors=True)
    os.makedirs("output_4_deal_inter")
    shutil.rmtree("output_5_inter_pep", ignore_errors=True)
    os.makedirs("output_5_inter_pep")
    shutil.rmtree("output_5_inter_cds", ignore_errors=True)
    os.makedirs("output_5_inter_cds")

    interpro_type = "Pfam"  # 改 - 选择不同的数据库（Pfam, CDD, Gene3D, MobiDBLite）
    pfam_num = {"DIR": ["PF03018"], "CCR": ["PF01370"], "COMT": ["PF00891"], "OMT": ["PF01596"], "P450": ["PF00067"]} # 改 - pfam 号的字典（可以是多个家族，每个家族对应一个列表）

    input_pep_files = glob.glob("input_pep_file/*")
    input_cds_files = glob.glob("input_cds_file/*")

    if os.path.exists("interproscan_pep_count.txt"):
        os.remove("interproscan_pep_count.txt")

    with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor:  # 改 - 线程池最大线程数
        for tsv_file in glob.glob("output_3_interproscan/*.tsv"):
            process_file = Deal_Inter(input_pep_files, input_cds_files, tsv_file, pfam_num, interpro_type)
            process_file.process_file()


    # 8 构建基因树文件，并画图
    gene_tree = GeneTree(glob.glob("output_5_inter_pep/*"))
    gene_tree.run()

    # 9 读取 block 文件，提取家族相关基因对
    gene_processor = GeneProcessor(glob.glob("output_6_seq_treeimg/*.fasta"), glob.glob(os.path.join("input_block_file", "*.rr.txt")))  # （可改） - block 文件夹
    gene_processor.run()

    # 10 去基因 pairs 进行去重，并合并文件
    gf_block_all = GFBlockAll(glob.glob(os.path.join("output_7_pair", "*.txt")))
    gf_block_all.run()

    # 11 判断基因对的 tandem < 50kb；
    gene_finder = GeneFinder(gff_files=glob.glob("input_gff_file/*.gff"), seq_treeimg_files=glob.glob("*_pair.txt"))
    gene_finder.run()

    end_time_all = time.time()
    print(f"End of run, total time is: {end_time_all - start_time_all} seconds")

