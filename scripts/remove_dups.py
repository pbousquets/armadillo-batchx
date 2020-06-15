#!/usr/bin/env python3
"""
This script filters duplicate mutations in Armadillo and discards those who are located in the flanking regions
of our regions of interest, as coverage issues may lead to FP.
"""

from sys import argv, stdin
flank_length = int(argv[1])

reads_dict = dict()

def rm_dups(line):
	chrom, pos, ID, REF, ALT, qual, filter, info, format, tumor, control = line.split()

	start = int(chrom.split(":")[1].split("-")[0])
	end = int(chrom.split(":")[1].split("-")[1])
	length = end - start
	
	if flank_length - 5 < int(pos) < length - flank_length + 5: #Remove the flanking extra regions. We leave 5bp because we still can trust those (not to far away) and belong to splicing sites.
		pass
	else: 
		return

	reads = set(read_name.split(","))
	length_reads = len(reads)
	if any(len(reads & other_mut_reads) > 0.5 * length_reads for other_mut_reads in reads_dict.values()): # if other mutation has 50% reads in common, consider it as duplicate
		return
	else:
		reads_dict[chrom + "_" + pos + "_" + REF + "_" + ALT] = reads
		print(line)

for line in stdin:
	line=line.strip()
	if line.startswith("#CHROM"):
		print('##Command=python3 %s' % (' '.join(argv)))
		print(line)
		continue
	elif line.startswith("##"):
		print(line)
	else:
		pass
	if len(line.split()) == 11:
		rm_dups(line)
	elif len(line.split()) == 22: 
		rm_dups("\t".join(line.split()[0:10]))
		rm_dups("\t".join(line.split()[10:]))
