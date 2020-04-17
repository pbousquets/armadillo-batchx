#!/usr/bin/env python3
"""
This script filters duplicate mutations in Armadillo and discards those who are located in the flanking regions
of our regions of interest, as coverage issues may lead to FP.
"""
from sys import argv, stdin
flank_length = int(argv[1])

reads_dict = dict()
header = ['#CHROM','POS','ID','REF','ALT', 'CHARACTERISTICS', "READS", "\n"]

for line in stdin:
	line=line.strip()
	if line.startswith("#"):
		print('##fileformat=VCFv4.2', '##Command=python3 %s' % (' '.join(argv)), line, '\t'.join(header), sep='\n', end='')
		continue

	chrom_TD, pos_TD, ID, REF_TD, ALT_TD, characteristics, read_name = line.split()
	start = int(chrom_TD.split(":")[1].split("-")[0])
	end = int(chrom_TD.split(":")[1].split("-")[1])
	length = end - start

	if flank_length - 5 < int(pos_TD) < length - flank_length + 5: #Remove the flanking extra regions. We leave 5bp because we still can trust those (not to far away) and belong to splicing sites.
		pass
	else: 
		continue

	key = chrom_TD+"_"+pos_TD+"_"+ID+"_"+REF_TD+"_"+ALT_TD
	reads = set(read_name.split(","))
	length_reads = len(reads)
	if any(len(reads & other_mut_reads) > 0.5 * length_reads for other_mut_reads in reads_dict.values()): # if other mutation has 50% reads in common, consider it as duplicate
		continue
	else:
		print(line)