import os
import glob

# 豆科
list_hh = ["Bva", "Adu", "Aed", "Aev", "Ahy", "Aip", "Amo", "Apr", "Car", "Cca", "Gma", "Dod", "Gso", "Lal", "Lan", "Lja", "Mal", "Mtr", "Psa", "Pvu", "Ssu", "Sto", "Tpr","Tsu", "Vra", "Vvi"]

# # 十字花科
# list = ["Asu","Ath","Aal","Ane","BjA","BjB","Bni","Bol","Bra" ,"Csa","Cru","Chi","Lma","Mpy","Aly","BnA","BnB","Rsa","Bvu","Tpa","Esa","Vvi"]
for line in list_hh:
	cmd = f"""
		mkdir {line}
		cp {line}*.png {line}
	"""
	os.system(cmd)




# files = glob.glob('*.fasta')

# for line in files:
# 	cmd = f"""/home/ubuntu/tools/kalign2/kalign -i {line} -o {line.split(".")[0]}_mul.fasta"""
# 	os.system(cmd)


# files_2 = glob.glob('*mul.fasta')
# for line in files_2:
# 	cmd = f"""iqtree2 -s {line} -m MFP -bb 1000 -nt AUTO"""
# 	os.system(cmd)
	

# from Bio import SeqIO


# for line in SeqIO.parse('gene.fasta', 'fasta'):
	
