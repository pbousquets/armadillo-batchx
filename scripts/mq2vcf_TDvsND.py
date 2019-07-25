#!/usr/bin/env python3

from sys import argv, stdin
from re import findall, match
from multiprocessing import Pool
from subprocess import check_output
import argparse

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description='Converts an MQ file into a VCF filtering by coverage')
	parser.add_argument(
	'-i', '--input', type=str, required=False, default='-', metavar='FILE',
	help='MQ file to convert')
	parser.add_argument(
	'-b', '--bam', type=str, required=True, metavar='FILE',
	help='Bam file to analyse')
	parser.add_argument(
	'-n', '--name', type=str, required=False, default='.',
	help='Value for the ID field. (default: %(default)s)')
	parser.add_argument(
	'-td', '--tumor_threshold', type=int, required=False, default=6, metavar='INT',
	help='Minimum number of mutated reads to consider a mutation in the tumor sample. (default: %(default)s)')
	parser.add_argument(
	'-nd', '--normal_threshold', type=int, required=False, default=2, metavar='INT',
	help='Maximum number of mutated reads to consider a mutation in the normal sample. (default: %(default)s)')
	parser.add_argument(
	'-f', '--full', action='store_true',
	help='Print all variants to follow how each is filtered in each step. No arguments required. (default: NULL)')
	parser.add_argument(
	'-gc', '--GCcutoff', type=int, required=False, default=80, metavar='INT',
	help='Max GC% per read. (default: %(default)s)')
	parser.add_argument(
	'-q', '--base_quality', type=int, required=False, default=30,
	help='Filter by base quality value (default: %(default)s)')
	parser.add_argument(
	'-r', '--reference_fasta', type=str, required=False, default='/data/databases/hg19_repetitive_exons_fasta/REF.fa',
	help='Reference fasta file (default: %(default)s)')
	parser.add_argument(
	'-c', '--normal_coverage', type=int, required=False, default=80,
	help='Minimal proportion (%) of normal coverage comparing to tumor coverage (default: %(default)s)')
	parser.add_argument(
	'-p', '--port', type=str, required=False, default='9001',
	help='Port for performing blat analysis (default: %(default)s)')
	parser.add_argument(
	'-t', '--threads', type=int, required=False, default=3,
	help='Number of threads to use (default: %(default)s)')
	return parser.parse_args()

def correct_variants_list(variants):
	variants, indel_list = list(variants), []
	while '+' in variants or '-' in variants or '^' in variants or '$' in variants: #To convert the indels in . to count them only once
		index_list = (j for j, k in enumerate(variants) if k == '+' or k == '-' or k == '^' or k == '$') #j is the index and k the element in variants, so if the element is + or -, it returns its index
		n = next(index_list) #To focus on the first element
		if variants[n] == '^': #^ represents the start of the read and the following character shows the mapQ of the read, thus we delete one of the symbols to take into account only once that mapQ value
			variants[n:n+2] = ''
			continue
		elif variants[n] == '$': #$ represents the end of a read segment, but the read mapQ value doesn't appear, so we have to remove the symbol
			variants[n:n+1] = ''
			continue
		else: #indel
			try:
				size = int(str(variants[n+1])+str(variants[n+2]))+2 #In case the indel had more than 9 changes
			except ValueError:
				size = int(variants[n+1])+1 #The size of the indel, which is the element following the - or + symbol
			indel_list.append(''.join(variants[n:n+size+1]))
			variants[n-1:n+size+1] = ';' #To convert the indel pattern into a . (e.g. .+4ACGT > ;), thus there is only one symbol for each mapQ value. Note that before the indel pattern there is a . or , which is the symbol that has the mapQ value, so we also have to replace it
			continue
	return(variants, indel_list)

def most_common_variant(single_variants_list, full_variants_list, threshold):
	most_count = []
	alt = []
	for variant in single_variants_list:
		count = full_variants_list.count(variant)
		if count >= threshold:
			most_count.append(str(count))
			alt.append(variant)
		else:
			pass
	return(alt)

def get_mut_reads(ref, alt, variants, reads):
	reads_ID_list = []
	if len(alt) == 1 and len(ref) == 1: #Substitutions
		look_for = alt
	else: #Indels
		look_for = ';'
	for base in zip(variants, reads):
		if base[0] == look_for: #To extract only those read names that contain the mutation we are looking for
			reads_ID_list.append(base[1])
		else:
			continue
	return(reads_ID_list)

def GC_filt(sequence):
	GC=int(sequence.count("G")) + int(sequence.count("C"))
	GC=GC/len(sequence)*100
	return GC

def samt_view(bamdir, coord, splitreads): # Samtools: get reads sequences
	samtools_command = ['samtools', 'view', bamdir, coord ]
	fastareads=''
	result = check_output(samtools_command).decode().splitlines()
	for readres in result:
		column=readres.split("\t")
		if column[0] in splitreads: #If the read is in the list of reads to analyze and isn't bad, analyze it
			GC=GC_filt(column[9]) #GC filter. If excesive GC content, remove the read
			if GC <= args.GCcutoff:
				fastareads=fastareads+">"+column[0]+"\n"+column[9]+"\n" #Store as fasta. "\n" allows blat to read it as multiline
	return fastareads

def blat_filter(blat_result): #Filter the blat result to remove the reads with 100% identity
	badreads=list()
	allreads=list()
	reads_left=[]
	for line in blat_result:
		column=line.split()
		try:
			ID=column[9]
			if ID not in allreads:
				allreads.append(ID) #store all IDs in a list
			if column[0]==column[10] and all(int(i) == 0 for i in column[1:8]) and int(column[11])==0 and column[12]==column[10] and ID not in badreads: #if perfect match, the change in the read is not a mutation..
				badreads.append(ID) #.. so append it to the bad reads list	
		except IndexError:
			pass
	for i in allreads:
		if i not in badreads:
			reads_left.append(i)
	nbadreads=len(badreads)
	return (reads_left, nbadreads)
	
def blat_search(fasta): # Blat search
	blat_command = ['gfClient', '-out=pslx', '-nohead','localhost', args.port, '', 'stdin', 'stdout']
	result = check_output(blat_command, input=fasta.encode()).decode().strip().split("\n")
	reads_left,nbadreads=blat_filter(result)
	return (reads_left, nbadreads)

def variantslist_correction(variants_list):
	variants = list(variants_list)
	while '+' in variants or '-' in variants or '^' in variants or '$' in variants: #To convert the indels in . to count them only once
		index_list = (j for j, k in enumerate(variants) if k == '+' or k == '-' or k == '^' or k == '$') #j is the index and k the element in variants, so if the element is + or -, it returns its index
		n = next(index_list) #To focus on the first element
		if variants[n] == '^': #^ represents the start of the read and the following character shows the mapQ of the read, thus we delete one of the symbols to take into account only once that mapQ value
			variants[n:n+2] = ''
			continue
		elif variants[n] == '$': #$ represents the end of a read segment, but the read mapQ value doesn't appear, so we have to remove the symbol
			variants[n:n+1] = ''
			continue
		else: #indel
			try:
				size = int(str(variants[n+1])+str(variants[n+2]))+2 #In case the indel had more than 9 changes
			except ValueError:
				size = int(variants[n+1])+1 #The size of the indel, which is the element following the - or + symbol
			variants[n-1:n+size+1] = str(variants[n+1]) #To convert the indel pattern into a . (e.g. .+4ACGT > 4), thus there is only one symbol for each mapQ value. Note that before the indel pattern there is a . or , which is the symbol that has the mapQ value, so we also have to replace it
			continue	
	return variants

def filter_changes(changesdic, allchanges, startdic, enddic, changesfrq):
	final_reads=list()
	for read in changesdic:
		errors=0
		bad=0
		for carrier in changesdic[read]: #Analyze the changes carried with the mutation
			position=int(carrier.split("_")[0])
			hits=int(allchanges.count(carrier)) #Frquency of the carrier in mut reads
			phasable_reads=0 #Reads that potentially could contain it (directly depends on its coordinates)
			for newread in changesdic: #Search how many reads could contain this mutation (i.e., the mutation is between its start and end)
				if position >= int(startdic[newread]) and position <= int(enddic[newread]):
					phasable_reads += 1
				else:
					pass
			if phasable_reads > 1 and hits == 1: #The change only appears once: probably an error 
				errors += 1
			elif 1 < hits < phasable_reads: #More than one read contains that change but not all: probably just an haplotype
				bad += 1
		if bad == 0 and errors < 4:
			final_reads.append(read)
	return final_reads

def parse_mq(pileup, reads_left, changesdic, changesfrq):
	startdic={}
	enddic={}
	allchanges=[]
	for newline in pileup:
		newline=newline.strip()
		column=newline.split()		
		chrom, pos, ref, coverage, variants, reads = column[0], column[1], column[2], column[3], column[4].upper(), column[6].split(',')
		for element in reads_left:
			if element in reads:
				if element not in startdic: #store the read start
					startdic[element]=pos
				else:
					pass
				enddic[element]=pos
				idx=reads.index(element) #store the read index to find its corresponding base in the variants list
				variants=variantslist_correction(variants)	
				if "," != variants[idx] and "." != variants[idx]: #Store everything but dots and commas (only keep the changes) for each read
					changesdic[element].append(pos+"_"+str(variants[idx]))
					allchanges.append(pos+"_"+str(variants[idx]))
					changesfrq[pos+"_"+str(variants[idx])]=variants.count(variants[idx])
				else:
					pass
			else:
				pass
	final_reads=filter_changes(changesdic, allchanges, startdic, enddic, changesfrq)
	return final_reads

def main_function(line):
	column = line.strip().split()
	chrom, pos, ref = column[0], column[1], column[2]
	coverage_td, variants_td, qual_td, reads_td = int(column[3]), column[4].upper(), column[5], column[6].split(',')
	coverage_nd, variants_nd = int(column[7]), column[8].upper()

	## Filter by minimun normal coverage
	if coverage_nd >= ((coverage_td*args.normal_coverage)/100):

		## Correct the variants lists and get the indels for the variants list
		corrected_variants_td, indel_list_td = correct_variants_list(variants_td)
		corrected_variants_nd, variants_list_nd = correct_variants_list(variants_nd)

		## Filter by base quality in the tumor sample
		new_corrected_variants_td, new_qual_td, new_reads_td, variants_list_td = [], [], [], []
		indel_count = 0
		for element in zip(corrected_variants_td, qual_td, reads_td):
			base_qual = ord(element[1])-33
			if base_qual >= args.base_quality:
				new_corrected_variants_td.append(element[0])
				new_qual_td.append(element[1])
				new_reads_td.append(element[2])
				if element[0] == ';':
					variants_list_td.append(indel_list_td[indel_count]) #To select only those indels that have a base quality greater than the cutoff
				else:
					pass
			else:
				pass
			if element[0] == ';':
				indel_count += 1
			else:
				pass
		## Get the possible mutations for this position
		variants_list_td += findall(r'(?<!\^)[ACGT]', ''.join(new_corrected_variants_td))
		variants_list_nd += findall(r'(?<!\^)[ACGT]', ''.join(corrected_variants_nd))

		## Count the frequency of variants
		variants_list_td_set = set(variants_list_td)
		if len(variants_list_td_set) == 0:
			pass
		else:
			alts_td = most_common_variant(variants_list_td_set, variants_list_td, args.tumor_threshold)
			alts_nd = most_common_variant(variants_list_td_set, variants_list_nd, args.normal_threshold+1)

			## Remove mutations that appear in normal sample
			real_alts_td = []
			for alt in alts_td:
				if alt not in alts_nd:
					real_alts_td.append(alt)
					continue
				else:
					if args.full:
						vsND=open("highND_lost.vcf", "a+")
						vsND.write(chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref+"\t"+alt+"\n")
						vsND.close()
					else:
						pass
			else:
				pass

			## Prepare the mutations to be printed
			ref_list = []
			if len(real_alts_td) == 0:
				pass
			else:
				for i in range(0, len(real_alts_td)):
					if len(real_alts_td[i]) == 1:
						ref_list.append(ref)
						continue
					else:
						nucleotides = ''.join(findall(r'[A-Z]', real_alts_td[i]))
						indel_deletion = match('-', real_alts_td[i])
						if indel_deletion is not None: #Deletions
							real_alts_td[i] = ref
							ref_list.append(ref + nucleotides)
							continue
						else: #Insertions
							real_alts_td[i] = ref + nucleotides
							ref_list.append(ref)
							continue
				for element in zip(ref_list, real_alts_td):
					real_reads_td = get_mut_reads(element[0], element[1], new_corrected_variants_td, new_reads_td)
					#Launch polyfilter
					bam=args.bam
					pos=int(pos)
					coord=chrom+":"+str(pos)+"-"+str(pos)
					reads_fasta=samt_view(bam, coord, real_reads_td)
					if len(reads_fasta)==0: #If no reads left, stop
						pass
					else:
						reads_left,nbadreads=blat_search(reads_fasta) #Launch blat to remove perfect reads
						if nbadreads < 4 and len(reads_left) >= args.tumor_threshold:
							changesdic={}
							changesfrq={}
							goodreads=[]
							if pos >= 200: #Phase the variants to check if all the reads containing the mutation come from the same region
								pileup_command=["samtools", "mpileup", "--output-QNAME", "-q", "0", "-Q", str(args.base_quality), "-R", "-f",  args.reference_fasta, bam, "-r", chrom+":"+str(int(pos)-200)+"-"+str(int(pos)+200)]
							else:
								pileup_command=["samtools", "mpileup", "--output-QNAME", "-q", "0", "-Q", str(args.base_quality), "-R", "-f",  args.reference_fasta, bam, "-r", chrom+":0-"+str(int(pos)+200)]

							pileup = check_output(pileup_command).decode().splitlines()
							for read in reads_left:
								changesdic[read]=[]
							final_reads=parse_mq(pileup, reads_left, changesdic, changesfrq)

							if len(final_reads) >= args.tumor_threshold:
								print(chrom, pos, args.name, element[0], element[1], ','.join(final_reads),"\n", sep='\t', end='')
							else:
								if args.full:
									f=open("phasing_lost.vcf", "a+")
									f.write(chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]+"\n")
									f.close()
						else:
							if args.full:
								f=open("blat_lost.vcf", "a+")
								f.write(chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]+"\n")
								f.close()
							else:
								pass
	else:
		f=open("total_coverage_lost.vcf", "a+")
		f.write(chrom+"\t"+pos+"\t"+str((100*coverage_nd/coverage_td))+"\n")
		f.close()

if __name__ == '__main__':
	args = parse_args()

	arg_list = []
	for arg in sorted(vars(args)):
		arg_list.append('--'+arg)
		arg_list.append(str(getattr(args, arg)))
	## VCF header
	header = ['#CHROM','POS','ID','REF','ALT','READS_NAME']
	print('##fileformat=VCFv4.2', '##Command=python3 %s %s' % (argv[0], ' '.join(arg_list)), '\t'.join(header), sep='\n')

	## Using multiprocessing.Pool each task is run in a differente thread with their own variables, and the result is given only when all the task are completed.
	pool = Pool(processes=args.threads)
	## Initialize the pool of threads
	if args.input == '-':
		pool.imap(func=main_function, iterable=stdin)
	else:
		pileup = open(args.input)
		pool.imap(func=main_function, iterable=pileup)
	pool.close()
	pool.join()

