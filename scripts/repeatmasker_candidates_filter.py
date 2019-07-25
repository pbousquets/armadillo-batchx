#!/usr/bin/env python3
#This file reads finalres. For each exon, it searches its coordinate (and its similar ones) in our repeatmaskerDB.
#A new column is added to finalres. If there are no repeats, a "-" is added, else the rep family is annotated. 
#Usage: cat unaan_candidates.vcf | grep -v ^#| python3 ${repeatmasker} ${repeatsDB}
from sys import argv, stdin
import tabix
repeatsDB=tabix.open(argv[1])
TDcutoff=int(argv[2])

#print header
header = ['#CHROM','POS','ID','REF','ALT','REPEATS', "READS"]
print('##fileformat=VCFv4.2', '##Command=python3 %s' % (' '.join(argv)), '\t'.join(header), sep='\n', end='\n')

readslist=[]
dic={}

def annotate(readslist, dic, line):
	chrom_TD, pos_TD, ID, REF_TD, ALT_TD, read_name = line.split() 
	bad=0
	good=list()
	reads=read_name.split(",")
	key=chrom_TD+"_"+pos_TD+"_"+ID+"_"+REF_TD+"_"+ALT_TD #We need a uniq key for each mutation
	##REMOVE DUPLICATES
	if key not in dic and read_name not in readslist: #Make a key for each mutation
		for element in reads:
			if element not in readslist:
				readslist.append(element) #By using that list, we avoid repetition of the same mutation (a uniq mutation may appear several times as they could appear in all repetitive exons)
				good.append(element)
			else:
				pass
		if len(good) >= TDcutoff:
			dic[key]=[",".join(good),list()] 
	##ANNOTATE REPEATS
			try:
				query=repeatsDB.querys(chrom_TD) #Do the query of the mut region 
				for element in query:
					starts=element[7].split(",")
					ends=element[8].split(",")
					for i in range(0, len(starts)): 
						if int(pos_TD) > int(starts[i]) and int(pos_TD) < int(ends[i]) and element[5] not in dic[key][1]: #If there's a repeat in there, annotate it.
							dic[key][1].append(element[5])
			except tabix.TabixError:
				pass
		else:
			pass
	else:
		pass
	return(readslist, dic)

for line in stdin:
	line=line.strip()
	if line.startswith("#"):
		pass
	else:
		if len(line.split())==6:
			readslist, dic=annotate(readslist, dic, line)
		elif len(line.split())==12:
			readslist, dic=annotate(readslist, dic, "\t".join(line.split()[0:6]))
			readslist, dic=annotate(readslist, dic, "\t".join(line.split()[6:]))
for i in dic:
	if len(dic[i][1])==0:
		dic[i][1]="-"
	else:
		pass
	print("\t".join(i.split("_")), dic[i][1], dic[i][0], sep="\t") #Print the dictionary
