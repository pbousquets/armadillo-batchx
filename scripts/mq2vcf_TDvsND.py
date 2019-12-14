#!/usr/bin/env python3
from sys import argv, stdin
from re import findall, match
from multiprocessing import Pool
from subprocess import check_output
import argparse
import statistics
import pandas as pd

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description = 'Converts an MQ file into a VCF filtering by coverage')
    parser.add_argument(
    '-c', '--normal_coverage', type = int, required = False, default = 80, metavar = 'INT',
    help = 'Minimal proportion (%) of normal coverage in regard to tumor coverage (default: %(default)s)')
    parser.add_argument(
    '-cb', '--control_bam', type = str, required = True, metavar = 'FILE',
    help = 'Tumor miniBAM to analyse')
    parser.add_argument(
    '-cov', '--genome_coverage', type = int, required = False, default = 30, metavar = 'INT',
    help = 'Tumor genome coverage')
    parser.add_argument(
    '-f', '--full', action = 'store_true', default = True,
    help = 'Print all variants to follow how each is filtered in each step. No arguments required. (default: NULL)')
    parser.add_argument(
    '-e', '--max_errors', type = int, required = False, default = 2, metavar = 'INT',
    help = 'Maximum context errors allowed in regard to read length(%%). (default: %(default)s %%)')
    parser.add_argument(
    '-se', '--sequencing_error', type = float, required = False, default = 0.0003, metavar = 'FLOAT',
    help = 'Sequencing error rate (default: %(default)s)')
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
        elif count_TD >= threshold_TD and count_ND > args.normal_max: #Store a list of those that were discarded due to high freq in control (likely SNPs)
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
        if column[0] in splitreads: #If the read is in the list of reads to analyse and isn't bad, analyse it
            GC = GC_filt(column[9]) #GC filter. If excesive GC content, remove the read
            if GC <= args.GCcutoff:
                fastareads = fastareads+">"+column[0]+"\n"+column[9]+"\n" #Store as fasta. "\n" allows blat to read it as multiline
    return fastareads

def blat_filter(blat_result): #Filter the blat result to remove the reads with 100% identity
    badreads = list()
    allreads = list()
    reads_left = list()
    for line in blat_result:
        column = line.strip().split("\t")
        try:
            ID = column[9]
            if ID not in allreads:
                allreads.append(ID) #store all IDs in a list
            if column[0] == column[10] and all(int(i) == 0 for i in column[1:8]) and int(column[11]) == 0 and column[12] == column[10] and ID not in badreads: #if perfect match, the change in the read is not a mutation..
                badreads.append(ID) #.. so append it to the bad reads list
            if column[1] > 1 and ID not in badreads: #if there are more than 1 changes, it unlikely will be a mutation
                badreads.append(ID) #.. so append it to the bad reads list
        except IndexError:
            pass
    reads_left = list(set(allreads)-set(badreads))
    nbadreads = len(badreads)
    return(reads_left, nbadreads)

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

def init_msa_matrix(reads_list, corrected_variants_list, pos):
    msa_data = dict()
    for i in range(0, len(reads_list)):
        msa_data[reads_list[i]] = corrected_variants_list[i]

    msa = pd.DataFrame.from_dict(msa_data, orient='index', columns = [pos])

    return(msa)

def msa_matrix(msa, reads_list, corrected_variants_list, pos):
    msa_data = dict()
    for i in range(0, len(reads_list)):
        msa_data[reads_list[i]] = corrected_variants_list[i]

    ### Append a new column to the matrix
    rownames = set(msa.index)
    to_complete = set(reads_list) - rownames #Add the rows for those read that started in the current position
    ended = rownames - set(reads_list) #Store those reads that already ended

    for each in to_complete: #Add new reads that didn't appear in previous positions. Add "*" in those positions for those reads.
        msa.loc[each] = ["*"] * msa.shape[1] #Add an asterisc in the previous positions of the reads that start in this position
        msa_data[each] = corrected_variants_list[reads_list.index(each)]

    for each in ended: #Add '*' for the ended ones
        msa_data[each] = "*"

    msa[pos] = pd.Series(msa_data)

    return(msa)

def pileup_to_msa(pileup):
    first = True
    noise = list()
    for line in pileup:
        column = line.strip().split()
        chrom, pos, ref = column[0], column[1], column[2]
        coverage_td, variants_td, reads_td = int(column[3]), column[4].upper(), column[6].split(',')
        coverage_nd, variants_nd, reads_nd = int(column[7]), column[8].upper(), column[10].split(',')

        ## Correct the variants lists and get the indels for the variants list ##
        corrected_variants_td, indel_list_td = correct_variants_list(variants_td)
        corrected_variants_nd, indel_list_nd = correct_variants_list(variants_nd)

        if first:
            td_msa = init_msa_matrix(reads_td, corrected_variants_td, pos)
            nd_msa = init_msa_matrix(reads_nd, corrected_variants_nd, pos)
            first = False
        else:
            td_msa = msa_matrix(td_msa, reads_td, corrected_variants_td, pos)
            nd_msa = msa_matrix(nd_msa, reads_nd, corrected_variants_nd, pos)

        noise.append(len([i for i, x in enumerate(corrected_variants_td) if x not in [",", "."]]))

    if "*" in td_msa.index: #Remove "*" as a row. It appears when in some positions there aren't reads and the pileup displays an asterisc.
        td_msa.drop("*", axis=0, inplace=True)

    if "*" in nd_msa.index: #Remove "*" as a row. It appears when in some positions there aren't reads and the pileup displays an asterisc.
        nd_msa.drop("*", axis=0, inplace=True)

    mean_noise = statistics.mean(noise)
    stdev_noise = statistics.stdev(noise)

    return(td_msa, nd_msa, mean_noise, stdev_noise)

def extract_context(msa, mut_pos, mut_base):
    context = set()
    context_frq = dict()
    bad_positions = 0
    context_error = 0
    mut_pos = str(mut_pos)
    only_mut_msa = msa[msa[mut_pos] == mut_base] #Filter the msa to keep only the mutant reads

    for pos in only_mut_msa.columns:
        if pos == mut_pos: #Don't consider the mutation as part of the context
            continue

        for variant in set(only_mut_msa[pos].tolist()):
            if variant in [",", ".", "*"]:
                continue
            variant_actual_frq = only_mut_msa[pos].tolist().count(variant) #"*" means the read doesn't map there, so we cannot count that one.
            context_discrepancies = len(only_mut_msa[pos].tolist()) - only_mut_msa[pos].tolist().count("*") - variant_actual_frq
            if variant_actual_frq == 1 and context_discrepancies > 0: #if the variant appears just once, it's just a sequencing error
                context_error += 1
                continue

            if context_discrepancies == 0:
                context.add(pos+"_"+variant)
                context_frq[pos+"_"+variant] = variant_actual_frq
            elif context_discrepancies == 1: #One discrepant read may be due just a sequencing error
                context.add(pos+"_"+variant)
                context_frq[pos+"_"+variant] = variant_actual_frq
                context_error += 1
            else:
                bad_positions += 1
    return(context, bad_positions, context_error, context_frq)

def analyse_tumor_context(msa, context, context_frq):
    unreliable_variants = 0
    for variant in context:
        ctxt_pos, ctxt_base = variant.split("_")
        context_variant_frq = int(context_frq[variant])
        total_variant_frq = msa[msa[ctxt_pos] == ctxt_base].shape[0] #Count the total context base frequency in the position

        if total_variant_frq <= context_variant_frq + 2: #If the context base only appears in the mutant reads, it's not realiable. Ask for more a couple of non-mutant reads
            unreliable_variants += 1
    return(unreliable_variants)

def analyse_control(msa, context, mut_pos, mut_base):
    mut_pos = str(mut_pos)
    ctxt_control_matches = dict()
    control_unreliable_variants = 0

    for variant in context:
        ctxt_pos, ctxt_base  = variant.split("_")
        ctxt_control_matches[variant] = msa[(msa[mut_pos] != mut_base) & (msa[mut_pos] != '*') & (msa[ctxt_pos] == ctxt_base)].shape[0]

    for variant in ctxt_control_matches.keys():
        control_unreliable_variants += 1 if ctxt_control_matches[variant] < 2 else control_unreliable_variants #If there are no reads with the context in the control, it's bad too

    return(control_unreliable_variants)

def analyse_context(chrom, mut_pos, mut_base):
    ## Extract the pileup ##
    pileup_start = mut_pos - args.read_length if mut_pos >= args.read_length else 1 #Start can't be negative
    pileup_end = mut_pos + args.read_length
    pileup_command = ["samtools", "mpileup", "--output-QNAME", "-q", "0", "-Q", str(args.base_quality), "-R", "-f",  args.reference_fasta, args.tumor_bam, args.control_bam, "-r", chrom+":"+str(pileup_start)+"-"+str(pileup_end)]
    pileup = check_output(pileup_command).decode().splitlines() #run samtools mpileup

    ## Convert the pileup into a msa ##
    td_msa, nd_msa, mean_noise, stdev_noise = pileup_to_msa(pileup)

    ## Analyse the tumor msa to extract the context ##
    context, bad_positions, context_error, context_frq = extract_context(td_msa, mut_pos, mut_base)

    if bad_positions > 0.2 * len(context): #We allow one discrepant position for every 5 context variants (20%). If greater, we'll discard de variant.
        ctxt_unreliable_variants = 99 #Set values for ctxt_unreliable_variants and control_ctxt_unreliable_variants so we can end the function without crushing.
        control_ctxt_unreliable_variants = 99
    else:
        ## Find the context in the normal sample and in unmutated tumoral reads ##
        ctxt_unreliable_variants = analyse_tumor_context(td_msa, context, context_frq)
        control_ctxt_unreliable_variants = analyse_control(nd_msa, context, mut_pos, mut_base)

    return(ctxt_unreliable_variants, control_ctxt_unreliable_variants, len(context), context_error, bad_positions, mean_noise, stdev_noise)


def main_function(line):
    column = line.strip().split()
    chrom, pos, ref = column[0], column[1], column[2]
    coverage_td, variants_td, qual_td, reads_td = int(column[3]), column[4].upper(), column[5], column[6].split(',')
    coverage_nd, variants_nd = int(column[7]), column[8].upper()
    ## Filter by minimun normal coverage
    if coverage_nd < ((coverage_td*args.normal_coverage)/100) or coverage_td < args.genome_coverage * 1.5: #If the genome is at 30x, at least we require 45x in each position as it should be repetitive
        if args.full:
            string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+"*"+"\t"+"*"
            if coverage_td < args.genome_coverage * 1.5:
                print_log(string, "low_tumor_coverage("+str(coverage_td)+"reads)")
            else:
                cov_ND = 100*coverage_nd/coverage_td
                print_log(string, "low_control_coverage("+str(int(cov_ND))+"%)")
        else:
            pass
        return
    else: #The mutation is good. Keep reading
        pass

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
        return # Stop if there are no mutations
    else:
        pass

    real_alts_td, bad_alts = most_common_variant(variants_list_td_set, variants_list_td, variants_list_nd, args.tumor_threshold, args.normal_max, coverage_td, coverage_nd)

    ## Annotate alts with frequency over threshold in control sample.
    if args.full:
        for alt in bad_alts:
            string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref+"\t"+alt
            print_log(string, "germline_change")
    else:
        pass

    ## Prepare the mutations to be printed
    if len(real_alts_td) == 0: #Stop if there are no mutations left
        return
    else:
        pass

    ref_list = list()

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
        pos = int(pos)
        coord = chrom+":"+str(pos)+"-"+str(pos)
        reads_fasta = samt_view(args.tumor_bam, coord, real_reads_td)
        if len(reads_fasta) == 0: #If no reads left, stop
            continue
        else:
            pass

        reads_left,nbadreads = blat_search(reads_fasta) #Launch blat to remove perfect reads
        if len(reads_left) < args.tumor_threshold: #If there are enough reads, go on
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "blat_perfect_match")
                continue
            else:
                continue
        else:
            pass

        ##Phase the variants to check if all the reads containing the mutation come from the same region
        mut_base = element[1]
        ctxt_unreliable_variants, control_ctxt_unreliable_variants, context_length, context_error, bad_positions, mean_noise, stdev_noise = analyse_context(chrom, pos, mut_base)

        if context_length != 0 and (ctxt_unreliable_variants/context_length > 0.2 or control_ctxt_unreliable_variants/context_length > 0.2):
            if args.full: #Remove the mutation if the context only exists in the mutant reads
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "context_errors")
                continue
            else:
                continue
        else:
            pass

        control_mutcov = variants_list_nd.count(mut_base)
        characteristics = [len(reads_left), coverage_td, nbadreads, control_mutcov, coverage_nd, mean_noise, stdev_noise]
        characteristics = list(map(str, characteristics))
        print(chrom, pos, args.name, element[0], element[1], ','.join(characteristics), ','.join(reads_left),"\n", sep = '\t', end = '')

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
