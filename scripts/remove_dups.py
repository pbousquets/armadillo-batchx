#!/usr/bin/env python3
"""
This script filters duplicate mutations in Armadillo and discards those who are located in the flanking regions
of our regions of interest, as coverage issues may lead to FP.
"""

from sys import argv, stdin
flank_length = int(argv[1])

reads_dict = dict()

def rm_dups(line):
	chrom, pos, ID, REF, ALT, qual, filter, info, format, tumor, control, read_name = line.split()

	reads = set(read_name.split(","))
	length_reads = len(reads)
	if any(len(reads & other_mut_reads) > 0.5 * length_reads for other_mut_reads in reads_dict.values()): # if other mutation has 50% reads in common, consider it as duplicate
		return
	else:
		reads_dict[chrom + "_" + pos + "_" + REF + "_" + ALT] = reads
		print(chrom, pos, ID, REF, ALT, qual, filter, info, format, tumor, control, sep = "\t")

for line in stdin:
	line=line.strip()
	if line.startswith("#CHROM"):
		print('##Command=python3 %s' % (' '.join(argv)))
		print(line)
		continue
	elif line.startswith("##"):
		print(line)
		continue
	else:
		pass
	
	if len(line.split()) == 12:
		rm_dups(line)
	elif len(line.split()) == 24: 
		rm_dups("\t".join(line.split()[0:11]))
		rm_dups("\t".join(line.split()[11:]))
