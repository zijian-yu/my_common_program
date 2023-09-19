#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Author :yuzijian1010@163.com
# @FileName :2.run_mc.py
# @Time :2023/5/31 18:41
# python3 2.run_mc.py
import re


# 1
def get_gene2pos(gff_file):
    """
    从gff文件中获取基因的位置信息
    :param gff_file: gff文件路径
    :return: 基因位置信息字典
    """
    gene2pos = {}
    with open(gff_file) as GG:
        for line in GG:
            arr = line.split()
            gene2pos[arr[5]] = arr[6]
    return gene2pos


def get_table(tc_gff_file, gene2pos):
    """
    从tc_gff文件中获取表格信息
    :param tc_gff_file: tc_gff文件路径
    :param gene2pos: 基因位置信息字典
    :return: 表格信息，tc基因编号字典，tc基因数量
    """
    table = []
    tcgenenum = 0
    tcgene2order = {}
    with open(tc_gff_file) as TG:
        for line in TG:
            arr = line.split()
            gene2pos[arr[5]] = arr[6]
            table.append([arr[5], " ", " "])
            tcgene2order[arr[5]] = tcgenenum
            tcgenenum += 1
    return table, tcgene2order, tcgenenum


def get_tc2Sp_ortho_regions(corr_file):
    """
    从corr文件中获取tc和Sp的正交区域信息
    :param corr_file: corr文件路径
    :return: tc和Sp的正交区域信息字符串
    """
    tc2Sp_ortho_regions = ""
    with open(corr_file) as CP:
        for line in CP:
            line = line.strip()
            tc2Sp_ortho_regions += line + "\t"
    return tc2Sp_ortho_regions


# 2
def get_Spgene2colnum(blockfile):
    """
    从block文件中获取Sp基因编号到列号的映射关系
    :param blockfile: block文件路径
    :return: Sp基因编号到列号的映射关系字典
    """
    with open(blockfile) as IN:
        colnum = 0
        Spgene2colnum = {}
        deal = ("+", " ", "\n", ">", "overlap")
        for line in IN:
            if line.startswith(deal):
                continue
            if line.startswith("t"):
                colnum = line.split()[4].strip()
                continue
            line = line.strip()
            id2, pos2, id1, pos1, orient = line.split()
            if id1 not in Spgene2colnum or Spgene2colnum[id1] < colnum:
                Spgene2colnum[id1] = colnum
    return Spgene2colnum


# 3
def get_blocklines(blockfile, gene2pos, tc2Sp_ortho_regions, Spgene2colnum, table, tcgene2order, type2col):
    """
    从block文件中获取块信息
    :param blockfile: block文件路径
    :param gene2pos: 基因位置信息字典
    :param tc2Sp_ortho_regions: tc和Sp的正交区域信息字符串
    :param Spgene2colnum: Sp基因编号到列号的映射关系字典
    :param table: 表格信息
    :param tcgene2order: tc基因编号字典
    :param type2col: 基因类型到列号的映射关系字典
    :return: 块信息列表
    """
    blocklines = []
    with open(blockfile) as IN:
        num = -1
        id1, id2, pos1, pos2, orient = None, None, None, None, None
        start1, start2, end1, end2 = None, None, None, None
        startid1, startid2, endid1, endid2 = None, None, None, None
        sp1, sp2, chr1, chr2 = None, None, None, None
        region2Spgene = {}
        colnum = 0
        deal = ("+", " ", "\n", ">", "overlap")
        for line in IN:
            if line.startswith(deal):
                continue
            if line.startswith("t"):
                colnum = line.split()[4].strip()
                colnum = colnum.replace("\n", "").replace("\r", "")
                num += 1
                continue
            if line.startswith(">"):
                pv = line.split()
                pv[3] = pv[3].replace("\n", "").replace("\r", "")
                end1 = gene2pos[id1]
                end2 = gene2pos[id2]
                endid1 = id1
                endid2 = id2
                if float(pv[3]) > 0.001 or endid2.startswith("Tc00") or colnum <= 4:
                    blocklines = []
                    num = -1
                    continue
            t2g = tc2Sp_ortho_regions.split()
            for i in range(len(t2g)):
                blocklines.append(line.strip())
                line = line.strip()
                id2, pos2, id1, pos1, orient = line.split()
                if len(blocklines) == 1:
                    sp1 = id1[:2]
                    sp2 = id2[:2]
                    chr1 = id1.split("g")[0]
                    chr2 = id2.split("g")[0]
                    start1 = gene2pos[id1]
                    start2 = gene2pos[id2]
                    startid1 = id1
                    startid2 = id2
                type, tcchr, tcseg_start, tcseg_end, Spchr, Spseg_start, Spseg_end = re.split(':|-', t2g[i])
                is2insert = 0
                for n in range(len(blocklines)):
                    id2, pos2, id1, pos1, orient = blocklines[n].split()
                    region = "_".join([Spchr, Spseg_start, Spseg_end, tcchr, tcseg_start, tcseg_end])
                    if region in region2Spgene:

                        Spgenes = region2Spgene["_".join([Spchr, Spseg_start, Spseg_end, tcchr, tcseg_start, tcseg_end])].split("\t")

                        for g in range(len(Spgenes)):
                            if id1 == Spgenes[g]:
                                is2insert -= 1
                    if Spgene2colnum.get(id1, 0) > colnum:
                        is2insert -= 1
                    if end1 is None:
                        end1 = 0
                    if isinstance(Spseg_start, str):
                        Spseg_start = int(Spseg_start)
                    if isinstance(Spseg_end, str):
                        Spseg_end = int(Spseg_end)
                    if isinstance(start1, str):
                        start1 = int(start1)
                    if isinstance(end1, str):
                        end1 = int(end1)
                    if isinstance(start2, str):
                        start2 = int(start2)
                    if isinstance(end2, str):
                        end2 = int(end2)
                    if (int(colnum) + int(is2insert)) >= 4 and tcchr == chr2 and Spchr == chr1 and int(min(start1, end1)) > int(Spseg_start)-10 and int(max(start1, end1)) < int(Spseg_end)+10:
                        if end2 is None:
                            end2 = 0
                        if isinstance(tcseg_start, str):
                            tcseg_start = int(tcseg_start)
                        if isinstance(tcseg_end, str):
                            tcseg_end = int(tcseg_end)
                        if (int(tcseg_start) + int(tcseg_end) > 0 and min(start2, end2) > int(tcseg_start)-10 and max(start2, end2) < int(tcseg_end)+10) or int(tcseg_start) + int(tcseg_end) == 0:
                            for n in range(len(blocklines)):
                                id2, pos2, id1, pos1, orient = blocklines[n].split()
                                table[tcgene2order[id2]][type2col[type]] = id1
                                if "_".join([Spchr, str(Spseg_start), str(Spseg_end), tcchr, str(tcseg_start), str(tcseg_end)]) not in region2Spgene:
                                    region2Spgene["_".join([Spchr, str(Spseg_start), str(Spseg_end), tcchr, str(tcseg_start), str(tcseg_end)])] = id1
                                else:
                                    region2Spgene["_".join([Spchr, str(Spseg_start), str(Spseg_end), tcchr, str(tcseg_start), str(tcseg_end)])] += "\t"+id1
                    blocklines = []
                    num = -1
    return table


# 4
def write_to_file(table, gene2pos, tcgenenum, species_name, output_file):
    """
    将table中的数据写入到文件中
    :param table: 二维列表，存储着基因的信息
    :param gene2pos: 字典，存储着基因的位置信息
    :param tcgenenum: int，tc基因的数量
    :param species_name: str，物种名称
    :param output_file: str，输出文件名
    :return: None
    """
    with open(output_file, "w") as TB:
        lastgenes = ["", "", "", "", "", ""]
        lastpos = [-1, -1, -1, -1, -1, -1]
        for i in range(tcgenenum):
            thisgenes = []
            for j in range(3):
                thisgenes.append(table[i][j])
                if j == 0:
                    continue
                if species_name not in thisgenes[j]:
                    continue
                if lastgenes[j] == "":
                    lastgenes[j] = thisgenes[j]
                    lastpos[j] = i
                elif lastgenes[j].startswith(species_name):
                    if lastgenes[j][:4] != thisgenes[j][:4]:
                        lastgenes[j] = thisgenes[j]
                        lastpos[j] = i
                    elif abs(int(gene2pos[lastgenes[j]]) - int(gene2pos[thisgenes[j]])) > 50:

                        lastgenes[j] = thisgenes[j]
                        lastpos[j] = i
                    else:
                        for p in range(lastpos[j]+1, i):
                            table[p][j] = "."
                        lastpos[j] = i
                        lastgenes[j] = thisgenes[j]
        for i in range(tcgenenum):
            depth = 0
            for j in range(3):
                if species_name in table[i][j] or table[i][j] == ".":
                    depth += 1
                TB.write(table[i][j] + ",")
            TB.write(str(depth) + "\n")


# 5
def convert_file(infile:str, outfile:str) -> None:
    """
    将infile中的数据转换成outfile中的数据
    :param infile: str，输入文件名
    :param outfile: str，输出文件名
    :return: None
    """
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        for line in f_in:
            array = line.strip().split(',')
            for i in range(4):
                if not array[i].isalnum():
                    array[i] = "."
            f_out.write('\t'.join(array[:4]) + '\n')


# 6
def process_files(file1:str, file2:str, output_file:str) -> None:
    """
    从file1中提取出所有的基因名字，存储到hash表中，然后从file2中提取出所有在hash表中的基因对应的行，写入到output_file中
    :param file1: str，文件名
    :param file2: str，文件名
    :param output_file: str，输出文件名
    :return: None
    """
    import re
    hash = {}
    with open(file1, "r") as f1:
        for line in f1:
            line = re.sub(r"\s+", ".", line.strip())
            line = re.sub(r",", "\t", line)
            line = re.sub(r"\t\d", "", line)
            arr = line.split("\t")
            for a in arr:
                if re.search(r"\d+g\d+", a):
                    b = arr[0] + "\t" + a
                    hash[b] = 1

    with open(file2, "r") as f2, open(output_file, "w") as out:
        for line in f2:
            array = line.strip().split("\t")
            c = array[0] + "\t" + array[1]
            if c in hash:
                out.write(line)


if __name__ == '__main__':
    # 1
    tc_gff_file = "Nda.new.gff"  # 改 参考物种的gff
    Sp_gff_file = "Aco.new.gff"  # 改 比对物种的gff
    corr_file = "Nda_Aco.corr.txt"  # 改 像素文件 -- 1.tran.posion.py 程序得到的结果
    gene2pos = get_gene2pos(Sp_gff_file)
    table, tcgene2order, tcgenenum = get_table(tc_gff_file, gene2pos)
    tc2Sp_ortho_regions = get_tc2Sp_ortho_regions(corr_file)
    # 2
    blockfile = "Nda_Aco.block.rr.txt"
    Spgene2colnum = get_Spgene2colnum(blockfile)
    # 3
    type2col = {"A1": 1, "A2": 2, "A3": 3, "A4": 4, "A5": 5, "A6": 6}  # 改 几比几的关系
    get_blocklines(blockfile, gene2pos, tc2Sp_ortho_regions, Spgene2colnum, table, tcgene2order, type2col)
    # 4
    mc_file = "CIR.Nda.Aco.mc.2.txt"  # 改 生成的mc文件命名（自己命名）
    write_to_file(table, gene2pos, tcgenenum, "Aco", mc_file)  # 改 物种名字
    # 5
    outfile = ".".join(mc_file.split(".")[:-1]) + "tab.dot.txt"
    convert_file(mc_file, outfile)
    # 6
    process_files(mc_file, "Nda_Aco.blast", "Nda_Aco.new.blast")  #  改 blast 文件，提取的blast文件


