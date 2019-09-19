#!/usr/bin/env python3
from sys import argv, stdin
import tabix
repeatsDB=tabix.open(argv[1])
TDcutoff=int(argv[2]) #Max duplicates allowed

#print header
header = ['#CHROM','POS','ID','REF', 'ALT', 'CHARACTERISTICS','REPEATS', "READS"]
print('##fileformat=VCFv4.2', '##Command=python3 %s' % (' '.join(argv)), '\t'.join(header), sep='\n', end='\n')

readslist=[]
dic={}
charact_dic={}

def annotate(readslist, dic, charact_dic, line):
	chrom_TD, pos_TD, ID, REF_TD, ALT_TD, CHARACTERISTICS, read_name = line.split()
	bad=list()
	reads=read_name.split(",")
	key=chrom_TD+"_"+pos_TD+"_"+ID+"_"+REF_TD+"_"+"_"+ALT_TD #We need a unique key for each mutation
	##REMOVE DUPLICATES
	if key not in dic and read_name not in readslist: #Make a key for each mutation
		for element in reads:
			if element not in readslist:
				readslist.append(element) #By using that list, we avoid duplicates of the same mutation (a unique mutation may appear several times as they may appear in any repeat of its real region)
				good.append(element)
			else:
				bad.append(element)
		if len(bad)/len(reads)* 100 <= TDcutoff:
			dic[key]=[",".join(good),list()]
			charact_dic=CHARACTERISTICS
	##ANNOTATE REPEATS
			try:
				query=repeatsDB.querys(chrom_TD) #Do the query of the mut region
				for result in query:
					starts=result[7].split(",")
					ends=result[8].split(",")
					for i in range(0, len(starts)):
						if int(pos_TD) > int(starts[i]) and int(pos_TD) < int(ends[i]) and result[5] not in dic[key][1]: #If there's a repeat in there, annotate it.
							dic[key][1].append(result[5])
			except tabix.TabixError:
				pass
		else:
			pass
	else:
		pass
	return(readslist, dic, charact_dic)

for line in stdin:
	line = line.strip()
	if line.startswith("#"):
		pass
	else:
		if len(line.split()) == 7:
			readslist, dic, charact_dic = annotate(readslist, dic, line, "\t".join(line.split()[0:7]))
		elif len(line.split()) == 14:  #For some reason, sometimes the previous script writes two mutations in the same line. If it happens, just split it.
			readslist, dic, charact_dic = annotate(readslist, dic, charact_dic, "\t".join(line.split()[0:7]))
			readslist, dic, charact_dic = annotate(readslist, dic, charact_dic, "\t".join(line.split()[7:]))
for key in dic:
	if len(dic[key][1])==0:
		dic[key][1] = "-"
	else:
		pass
	print("\t".join(key.split("_")), charact_dic[key], [key][1], dic[key][0], sep="\t") #Print the dictionary
