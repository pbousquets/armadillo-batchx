#!/usr/bin/env python3

from sys import argv, stdin
from re import findall, match
from multiprocessing import Pool
from subprocess import check_output
import argparse
import statistics

def parse_args():
	"""Parse the input arguments, use '-h' for help"""
	parser = argparse.ArgumentParser(description = 'Converts an MQ file into a VCF filtering by coverage')
	parser.add_argument(
	'-c', '--normal_coverage', type = int, required = False, default = 80,
	help = 'Minimal proportion (%) of normal coverage in regard to tumor coverage (default: %(default)s)')
	parser.add_argument(
	'-cb', '--control_bam', type = str, required = True, metavar = 'FILE',
	help = 'Tumor miniBAM to analyse')
	parser.add_argument(
	'-f', '--full', action = 'store_true', default = True,
	help = 'Print all variants to follow how each is filtered in each step. No arguments required. (default: NULL)')
	parser.add_argument(
	'-e', '--max_errors', type = int, required = False, default = 2, metavar = 'INT',
	help = 'Maximum context errors allowed in regard to read length(%%). (default: %(default)s %%)')
	parser.add_argument(
	'-i', '--input', type = str, required = False, default = '-', metavar = 'FILE',
	help = 'MQ file to convert')
	parser.add_argument(
	'-gc', '--GCcutoff', type = int, required = False, default = 80, metavar = 'INT',
	help = 'Max GC% per read. (default: %(default)s)')
	parser.add_argument(
	'-n', '--name', type = str, required = False, default = '.',
	help = 'Value for the ID field. (default: %(default)s)')
	parser.add_argument(
	'-nt', '--normal_contamination', type = int, required = False, default = 15, metavar = 'INT',
	help = 'Contamination %% in normal sample. (default: %(default)s)')
	parser.add_argument(
	'-nm', '--normal_max', type = int, required = False, default = 3, metavar = 'INT',
	help = 'Maximum mutant reads allowed in the control. (default: %(default)s)')
	parser.add_argument(
	'-p', '--port', type = str, required = False, default = '9001',
	help = 'Port for performing blat analysis (default: %(default)s)')
	parser.add_argument(
	'-q', '--base_quality', type = int, required = False, default = 30,
	help = 'Filter by base quality value (default: %(default)s)')
	parser.add_argument(
	'-r', '--reference_fasta', type = str, required = True,
	help = 'Reference fasta file')
	parser.add_argument(
	'-rl', '--read_length', type = int, required = False, default = 151,
	help = 'Average read length (default: %(default)s)')
	parser.add_argument(
	'-t', '--threads', type = int, required = False, default = 3,
	help = 'Number of threads to use (default: %(default)s)')
	parser.add_argument(
	'-tb', '--tumor_bam', type = str, required = True, metavar = 'FILE',
	help = 'Tumor miniBAM to analyse')
	parser.add_argument(
	'-tt', '--tumor_threshold', type = int, required = False, default = 4, metavar = 'INT',
	help = 'Minimum number of mutated reads to consider a mutation in the tumor sample. (default: %(default)s)')
	return parser.parse_args()

def print_log(line, reason):
	f = open("discarded_variants.log", "a+")
	f.write(line+"\t"+reason+"\n")
	f.close()

def correct_variants_list(variants):
	variants, indel_list = list(variants), list()
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

def most_common_variant(single_variants_list, full_variants_list_TD, full_variants_list_ND, threshold_TD, normal_max, coverage_td, coverage_nd):
	alts = list()
	badalts = list()
	for variant in single_variants_list:
		count_TD = full_variants_list_TD.count(variant)
		count_ND = full_variants_list_ND.count(variant)
		if count_TD >= threshold_TD and count_ND <= args.normal_max : #We don't allow more than 3 reads in a mutation
			alts.append(variant)
		elif count_TD >= threshold_TD and count_ND > realND_threshold: #Store a list of those that were discarded due to high freq in control (likely SNPs)
			badalts.append(variant)
	return(alts, badalts)

def get_mut_reads(ref, alt, variants, reads):
	reads_ID_list = list()
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
	GC = int(sequence.count("G")) + int(sequence.count("C"))
	GC = GC/len(sequence)*100
	return GC

def samt_view(bamdir, coord, splitreads): # Samtools: get reads' sequences
	samtools_command = ['samtools', 'view', bamdir, coord ]
	fastareads = ''
	result = check_output(samtools_command).decode().splitlines()
	for readres in result:
		column = readres.split("\t")
		if column[0] in splitreads: #If the read is in the list of reads to analyze and isn't bad, analyze it
			GC = GC_filt(column[9]) #GC filter. If excesive GC content, remove the read
			if GC <= args.GCcutoff:
				fastareads = fastareads+">"+column[0]+"\n"+column[9]+"\n" #Store as fasta. "\n" allows blat to read it as multiline
	return fastareads

def blat_filter(blat_result): #Filter the blat result to remove the reads with 100% identity
	badreads = list()
	allreads = list()
	reads_left = list()
	for line in blat_result:
		column = line.split()
		try:
			ID = column[9]
			if ID not in allreads:
				allreads.append(ID) #store all IDs in a list
			if column[0] == column[10] and all(int(i) == 0 for i in column[1:8]) and int(column[11]) == 0 and column[12] == column[10] and ID not in badreads: #if perfect match, the change in the read is not a mutation..
				badreads.append(ID) #.. so append it to the bad reads list
		except IndexError:
			pass
	reads_left = list(set(allreads)-set(badreads))
	nbadreads = len(badreads)
	return (reads_left, nbadreads)

def blat_search(fasta): # Blat search
	blat_command = ['gfClient', '-out=pslx', '-nohead','localhost', args.port, '', 'stdin', 'stdout']
	result = check_output(blat_command, input = fasta.encode()).decode().strip().split("\n")
	reads_left,nbadreads = blat_filter(result)
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

def variants_extract(variants, reads, pos, start_dict, end_dict, variants_dict, mut_pos):
	variants = variantslist_correction(variants)
	indices = [i for i, x in enumerate(variants) if x not in [",", "."]]
	mut_prefix = str(mut_pos)+"_"
	for idx in indices:
		if reads[idx].startswith(mut_prefix):
			continue
		else:
			start_dict[reads[idx]] = pos if reads[idx] not in start_dict else start_dict[reads[idx]] #Store the start of any read
			end_dict[reads[idx]] = pos #Store the end of any read (the value will be replaced over and over until the last base of the read)
			if reads[idx] not in variants_dict:
				variants_dict[reads[idx]] = [str(pos)+"_"+str(variants[idx])]
			else:
				variants_dict[reads[idx]].append(str(pos)+"_"+str(variants[idx])) #Save the variants of each read in a dict
	return(start_dict, end_dict, variants_dict)

def context_extract(changes_dict, start_dict, end_dict):
	read_max_errors = args.max_errors*args.read_length/100
	read_max_errors = int(read_max_errors) + 1 if read_max_errors % 1 != 0 else read_max_errors #Round it up
	errors = 0
	context = set()
	different_context = set()
	changes_set = set()
	changes_list = list()

	for each_list in changes_dict.values(): # Merge all changes found in the reads into one set and create a list of all changes
		changes_set = changes_set | set(each_list)
		changes_list.extend(each_list)

	for change in changes_set:
		change_pos = int(change.split("_")[0])

		possible_reads = set([read for read, start in start_dict.items() if int(start) <= change_pos]) & set([read for read, end in end_dict.items() if int(end) <= change_pos]) #Phasable if the change is between start and end
		expected_change_frequency = len(possible_reads) #All reads that are can own the change should carry it.
		actual_change_frequency = changes_list.count(change)
		if expected_change_frequency - actual_change_frequency == 0: #All reads contain the change, so it's part of the context
			context.add(change)
		elif actual_change_frequency == 1 and expected_change_frequency > 1: #Only appears once --> error. Not important
			errors += 1
		elif expected_change_frequency - actual_change_frequency == 1: #One of them doesn't contain it. Probably an error but it belongs to the context
			errors += 1
			context.add(change)
		else: #Sometimes it appears, sometimes it doesn't --> different origin.
			different_context.add(change)

	#Maybe not all context variants appear in mutant reads just because of their start and end. Take it into account:
	reads_possible_context = {read:[variant for variant in context if int(variant.split("_")[0]) <= end_dict[read] and int(variant.split("_")[0]) >= start_dict[read]] for read in changes_dict.keys()} #The reads may not contain all context variants because of the start and end positions, so we take it into account

	##Define what reads are bad##
	errors_set = changes_set - context - different_context #Errors: variants that don't define the mutation context nor any other one.
	bad_reads = set(read for read, read_changes in changes_dict.items() if len(set(read_changes) & different_context) > 0) #Those that have a change that doesn't define the context of the mutation
	bad_reads = bad_reads & set(read for read, read_changes in changes_dict.items() if len(set(read_changes) & errors_set) > read_max_errors) #If errors are above threshold, remove the read as well
	bad_reads = bad_reads & set(read for read, read_changes in reads_possible_context.items() if len(set(read_changes) - set(changes_dict[read])) > 1 ) #Those that miss more than one element of context are bad too.

	final_reads = list(set(changes_dict.keys())-bad_reads)

	return(final_reads, bad_reads, context)

def context_analysis(pileup, reads_left, mut_pos):
	tumor_variants, tumor_startdict, tumor_enddict = (dict() for i in range(3))
	control_variants, control_startdict, control_enddict = (dict() for i in range(3))
	noise = list()
	for newline in pileup: # READ THE PILEUP
		column = newline.strip().split()
		chrom, pos, ref, coverage, tumor_variants_list, reads, control_coverage, control_variants_list, control_reads = column[0], column[1], column[2], column[3], column[4].upper(), column[6].split(','), column[7], column[8].upper(), column[10].split(',')
		noise.append(len([i for i, x in enumerate(tumor_variants_list) if x not in [",", "."]])) #Record how many variants there are in this position
		if control_variants_list in ['*'] or tumor_variants_list in ['*']: #Ignore that position if there's no coverage
			continue
		else:
			pass
		tumor_startdict, tumor_enddict, tumor_variants = variants_extract(tumor_variants_list, reads, pos, tumor_startdict, tumor_enddict, tumor_variants, mut_pos)
		control_startdict, control_enddict, control_variants = variants_extract(control_variants_list, control_reads, pos, control_startdict, control_enddict, control_variants, mut_pos)

	## GET THE CONTEXT ##
	mutantreads_changes_dict = dict((read, tumor_variants[read]) for read in reads_left) #Save the mutant reads left in a dictionary
	final_reads, badreads, context = context_extract(mutantreads_changes_dict, tumor_startdict, tumor_enddict) # Extracts the mutant reads with the same context

	## FIND THEORETICAL NON-MUTANT READS WITH SAME CONTEXT ##
	if len(context) == 0: #If there are no SNPs associated to the mutation:
		context_tumor_reads = len({read:[variant for variant in context if int(variant.split("_")[0]) <= tumor_enddict[read] and int(variant.split("_")[0]) >= tumor_startdict[read]] for read in tumor_variants.keys()})
		context_control_reads = len({read:[variant for variant in context if int(variant.split("_")[0]) <= control_enddict[read] and int(variant.split("_")[0]) >= control_startdict[read]] for read in control_variants.keys()})
	else:
		tumor_ctxtpossible_reads = {read:[variant for variant in context if int(variant.split("_")[0]) <= tumor_enddict[read] and int(variant.split("_")[0]) >= tumor_startdict[read]] for read in tumor_variants.keys()} #The reads may not contain all context variants because of the start and end positions, so we take it into account. variant.split("_")[0] is the change position
		context_tumor_reads=len([read for read, read_changes in tumor_variants.items() if set(tumor_ctxtpossible_reads.keys()).issubset(set(read_changes))])

		control_ctxtpossible_reads = {read:[variant for variant in context if int(variant.split("_")[0]) <= control_enddict[read] and int(variant.split("_")[0]) >= control_startdict[read]] for read in control_variants.keys()} #The reads may not contain all context variants because of the start and end positions, so we take it into account
		context_control_reads=len([read for read, read_changes in control_variants.items() if set(control_ctxtpossible_reads.keys()).issubset(set(read_changes))])

	stdev_noise = statistics.stdev(noise)
	mean_noise = statistics.mean(noise)

	return(final_reads, context, mean_noise, stdev_noise, context_tumor_reads, context_control_reads)

def main_function(line):
	column = line.strip().split()
	chrom, pos, ref = column[0], column[1], column[2]
	coverage_td, variants_td, qual_td, reads_td = int(column[3]), column[4].upper(), column[5], column[6].split(',')
	coverage_nd, variants_nd = int(column[7]), column[8].upper()
	## Store all variants in the region for context analysis if needed --> determinar si hace falta esto
	## Filter by minimun normal coverage
	if coverage_nd >= ((coverage_td*args.normal_coverage)/100):
		## Correct the variants lists and get the indels for the variants list
		corrected_variants_td, indel_list_td = correct_variants_list(variants_td)
		corrected_variants_nd, variants_list_nd = correct_variants_list(variants_nd)

		## Filter by base quality in the tumor sample
		new_corrected_variants_td, new_reads_td, variants_list_td = [], [], []
		indel_count = 0
		for element in zip(corrected_variants_td, qual_td, reads_td):
			base_qual = ord(element[1])-33
			if base_qual >= args.base_quality:  #If quality below threshold, this may be just a sequencing error
				new_corrected_variants_td.append(element[0])
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

			real_alts_td, bad_alts = most_common_variant(variants_list_td_set, variants_list_td, variants_list_nd, args.tumor_threshold, args.normal_max, coverage_td, coverage_nd)

			## Annotate alts with frequency over threshold in control sample.
			if args.full:
				for alt in bad_alts:
					string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref+"\t"+alt
					print_log(string, "germline_change")
			else:
				pass

			## Prepare the mutations to be printed
			ref_list = list()
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

					##Launch polyfilter
					tbam = args.tumor_bam
					cbam = args.control_bam
					pos = int(pos)
					coord = chrom+":"+str(pos)+"-"+str(pos)
					reads_fasta = samt_view(tbam, coord, real_reads_td)
					if len(reads_fasta) == 0: #If no reads left, stop
						pass
					else:
						reads_left,nbadreads = blat_search(reads_fasta) #Launch blat to remove perfect reads
						if len(reads_left) >= args.tumor_threshold: #If there are enough reads, go on
							##Phase the variants to check if all the reads containing the mutation come from the same region
							pileup_start = pos-args.read_length if pos >= args.read_length else 1 #Start can't be negative
							pileup_command = ["samtools", "mpileup", "--output-QNAME", "-q", "0", "-Q", str(args.base_quality), "-R", "-f",  args.reference_fasta, tbam, cbam, "-r", chrom+":"+str(pileup_start)+"-"+str(pos+args.read_length)]
							pileup = check_output(pileup_command).decode().splitlines() #run samtools mpileup

							final_reads, context, mean_noise, stdev_noise, context_tumor_reads, context_control_reads = context_analysis(pileup, reads_left, pos) #Analyse the context of each variant
							nbadreads += len(reads_left) - len(final_reads) #nbadreads equals those discarded with blat plus the ones discarded with context analysis

							if args.full and (context_control_reads < 2 or context_tumor_reads - len(final_reads) < 2): #Remove the mutation if the context only exists in the mutant reads
								string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
								print_log(string, "no_context_without_mut")
								continue
							else:
								pass

							corrected_control_reads_threshold = context_control_reads * len(final_reads)/context_tumor_reads * args.normal_threshold/100 #Use context reads as total coverage (it's more close to the real coverage than using the raw depth) to compute the MAF
							control_mutcov = variants_list_nd.count(element[1]) #Count frequency of mut in the control
							if corrected_control_threshold < control_mutcov and args.full: #Check if the mut coverage is over the threshold
								string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
								print_log(string, "too_many_mutreads_in_control("+str(control_mutcov)+")")
							elif corrected_control_threshold >= control_mutcov:
								pass #Go on, the mutation is good.
							else:
								continue # Too many mutant reads in the control, but don't print it because args.full is FALSE

							if len(final_reads) >= args.tumor_threshold:
								characteristics = [len(final_reads), coverage_td, nbadreads, context_tumor_reads, control_mutcov, coverage_nd, context_control_reads, mean_noise, stdev_noise]
								characteristics = list(map(str, characteristics))
								print(chrom, pos, args.name, element[0], element[1], ','.join(characteristics), ','.join(final_reads),"\n", sep = '\t', end = '')
							elif args.full and len(final_reads) < args.tumor_threshold:
								string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
								print_log(string, "many_bad_reads")
							else:
								pass
						elif args.full:
							string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
							print_log(string, "blat_perfect_match")
						else:
							pass
	elif args.full:
		string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+"*"+"\t"+"*"
		cov_ND = 100*coverage_nd/coverage_td
		print_log(string, "low_control_coverage("+str(cov_ND)+"%)")
	else:
		pass

if __name__ == '__main__':
	args = parse_args()

	arg_list = list()
	for arg in sorted(vars(args)):
		arg_list.append('--'+arg)
		arg_list.append(str(getattr(args, arg)))
	## VCF header
	header = ['#CHROM','POS','ID','REF','ALT', 'CHARACTERISTICS', 'READS_NAME']
	print('##fileformat = VCFv4.2', '##Command = python3 %s %s' % (argv[0], ' '.join(arg_list)), '\t'.join(header), sep = '\n')

	## Using multiprocessing.Pool each task is run in a differente thread with their own variables, and the result is given only when all the task are completed.
	pool = Pool(processes = args.threads)
	## Initialize the pool of threads
	if args.input == '-':
		pool.imap(func = main_function, iterable = stdin)
	else:
		pileup = open(args.input)
		pool.imap(func = main_function, iterable = pileup)
	pool.close()
	pool.join()
