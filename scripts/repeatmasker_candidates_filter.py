#This file reads finalres. For each exon, it searches its coordinate (and its similar ones) in our repeatmaskerDB.
#A new column is added to finalres. If there are no repeats, a "-" is added, else the rep family is annotated.
#Usage: cat unaan_candidates.vcf | grep -v ^#| python3 ${repeatmasker} ${repeatsDB}
from sys import argv, stdin
import tabix
repeatsDB=tabix.open(argv[1])
TDcutoff=int(argv[2])
flank_length=int(argv[3])
#print header
header = ['#CHROM','POS','ID','REF','ALT', 'CHARACTERISTICS', 'REPEATS', "READS"]
print('##fileformat=VCFv4.2', '##Command=python3 %s' % (' '.join(argv)), '\t'.join(header), sep='\n', end='\n')

readslist=[]
dic={}

def annotate(readslist, dic, line):
		##REMOVE DUPLICATES
		chrom_TD, pos_TD, ID, REF_TD, ALT_TD, read_name = line.strip().split()
		start=int(chrom_TD.split(":")[1].split("-")[0])
		end=int(chrom_TD.split(":")[1].split("-")[1])
		length=end-start
		if flank_length< int(pos_TD) < length - flank_length: #Remove the flanking extra regions. We leave 5bp because we still can trust those (not to far away) and belong to splicing sites.
				reads=read_name.split(",")
				bad=0
				key=chrom_TD+"_"+pos_TD+"_"+ID+"_"+REF_TD+"_"+ALT_TD #We need a uniq key for each mutation
				pos=int(chrom_TD.split(":")[1].split("-")[0])+ int(pos_TD)
				if key not in dic and read_name not in readslist: #Make a key for each mutation
						for element in reads:
								if element not in readslist:
										readslist.append(element) #By using that list, we avoid repetition of the same mutation (as a uniq mutation may appear several times as they should appear in all repetitive exons)
								else:
										bad += 1
						if bad < TDcutoff:
								dic[key]=[read_name,set()]
								##ANNOTATE REPEATS
								try:
										query=repeatsDB.querys(chrom_TD) #Do the query of the mut region
										for element in query:
												dic[key][1].add(element[3])
								except tabix.TabixError:
										pass
						else:
								if len(argv)==5:
										f=open("duplicates.vcf", "a+")
										f.write(line+"\n")
		else:
				if len(argv)==5:
						f=open("duplicates.vcf", "a+")
						f.write(line.split()+"\n")
		return(readslist, dic)

for line in stdin:
		line=line.strip()
		if line.startswith("#"):
				pass
		else:
				if len(line.split())==6:
						readslist, dic=annotate(readslist, dic, line)
				elif len(line.split())==12:
						readslist, dic=annotate(readslist, dic, "\t".join(line.split()[0:7]))
						readslist, dic=annotate(readslist, dic, "\t".join(line.split()[7:]))
for i in dic:
		if len(dic[i][1])==0:
				dic[i][1]="-"
		else:
				pass
		print("\t".join(i.split("_")), ",".join(dic[i][1]), dic[i][0], sep="\t") #Print the dictionary
