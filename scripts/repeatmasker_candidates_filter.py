#This file reads finalres. For each exon, it searches its coordinate (and its similar ones) in our repeatmaskerDB.
#A new column is added to finalres. If there are no repeats, a "-" is added, else the rep family is annotated.
#Usage: cat unaan_candidates.vcf | grep -v ^#| python3 ${repeatmasker} ${repeatsDB}
from sys import argv, stdin
import tabix
repeatsDB=tabix.open(argv[1])
TDcutoff=int(argv[2]) #Max duplicates allowed (%)
flank_length=int(argv[3])

readslist=[]
dic={}
charact_dic={}

def annotate(readslist, dic, line):
		##REMOVE DUPLICATES FIRST
		chrom_TD, pos_TD, ID, REF_TD, ALT_TD, characteristics, read_name = line.strip().split()
		start=int(chrom_TD.split(":")[1].split("-")[0])
		end=int(chrom_TD.split(":")[1].split("-")[1])
		length=end-start
		bad=list()

		if flank_length - 5 < int(pos_TD) < length - flank_length + 5: #Remove the flanking extra regions. We leave 5bp because we still can trust those (not to far away) and belong to splicing sites.
				reads=read_name.split(",")
				key=chrom_TD+"_"+pos_TD+"_"+ID+"_"+REF_TD+"_"+ALT_TD #We need a unique key for each mutation
				pos=int(chrom_TD.split(":")[1].split("-")[0])+ int(pos_TD)
				if key not in dic and read_name not in readslist: #Make a key for each mutation
						for result in reads:
								if result not in readslist:
										readslist.append(result) #By using that list, we avoid duplicates of the same mutation (a unique mutation may appear several times as they may appear in any repeat of its real region)
								else:
										bad.append(result)
						if len(bad)/len(reads)* 100 <= TDcutoff:
								dic[key]=[",".join(reads),list()]
								charact_dic[key]=characteristics
								##ANNOTATE REPEATS
								try:
										query=repeatsDB.querys(chrom_TD) #Do the query of the mut region
										for element in query:
												dic[key][1]=(element[3])
								except tabix.TabixError:
										pass
						else:
							pass
		else:
				pass
		return(readslist, dic, charact_dic)

#print header
header = ['#CHROM','POS','ID','REF','ALT', 'CHARACTERISTICS', 'REPEATS', "READS"]
print('##fileformat=VCFv4.2', '##Command=python3 %s' % (' '.join(argv)), '\t'.join(header), sep='\n', end='\n')

for line in stdin:
		line=line.strip()
		if line.startswith("#"):
				pass
		else:
				if len(line.split())==7:
						readslist, dic, charact_dic=annotate(readslist, dic, line)
				elif len(line.split())==14:
						readslist, dic, charact_dic=annotate(readslist, dic, "\t".join(line.split()[0:7]))
						readslist, dic, charact_dic=annotate(readslist, dic, "\t".join(line.split()[7:]))
for key in dic:
		if len(dic[key][1])==0:
				dic[key][1]="-"
		else:
				pass
		print("\t".join(key.split("_")), charact_dic[key], dic[key][1], dic[key][0], sep="\t") #Print the dictionary
