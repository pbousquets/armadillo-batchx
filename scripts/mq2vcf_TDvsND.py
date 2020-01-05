#!/usr/bin/env python3
from sys import argv, stdin
from re import findall, match
from multiprocessing import Pool
from subprocess import check_output
import argparse
import statistics
import pandas as pd
import numpy as np

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description = 'Converts an MQ file into a VCF filtering by coverage')
    parser.add_argument(
    '-cc', '--control_coverage', type = int, required = False, default = 80, metavar = 'INT',
    help = 'Minimal proportion (%) of control coverage in regard to tumor coverage (default: %(default)s)')
    parser.add_argument(
    '-cb', '--control_bam', type = str, required = True, metavar = 'FILE',
    help = 'Tumor miniBAM to analyse')
    parser.add_argument(
    '-tc', '--tumor_coverage', type = int, required = False, default = 30, metavar = 'INT',
    help = 'Tumor genome coverage')
    parser.add_argument(
    '-f', '--full', action = 'store_true', default = True,
    help = 'Print all variants to follow how each is filtered in each step. No arguments required. (default: NULL)')
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
    '-cm', '--control_max', type = int, required = False, default = 3, metavar = 'INT',
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

def most_common_variant(single_variants_list, full_variants_list_TD, full_variants_list_ND, threshold_TD, control_max, coverage_td, coverage_nd):
    alts = list()
    badalts = list()
    for variant in single_variants_list:
        count_TD = full_variants_list_TD.count(variant)
        count_ND = full_variants_list_ND.count(variant)
        if args.tumor_coverage > count_TD >= threshold_TD and count_ND <= args.control_max : #We don't allow more than 3 reads in a mutation
            alts.append(variant)
        elif count_TD >= threshold_TD and count_ND > args.control_max: #Store a list of those that were discarded due to high freq in control (likely SNPs)
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
    badreads, goodreads, allreads, reads_left = list(), list(), list(), list()
    for line in blat_result:
        column = line.strip().split("\t")
        try:
            ID = column[9]
            if ID not in allreads:
                allreads.append(ID) #store all IDs in a list
            if column[0] == column[10] and all(int(i) == 0 for i in column[1:8]) and int(column[11]) == 0 and column[12] == column[10] and ID not in badreads: #if perfect match, the change in the read is not a mutation..
                badreads.append(ID) #.. so append it to the bad reads list
                continue
        except IndexError:
            pass
    reads_left = list(set(allreads) - set(badreads))
    return(reads_left, len(badreads))

def blat_search(fasta): # Blat search
    blat_command = ['gfClient', '-out=pslx', '-nohead','localhost', args.port, '', 'stdin', 'stdout']
    result = check_output(blat_command, input = fasta.encode()).decode().strip().split("\n")
    reads_left,nbadreads = blat_filter(result)
    return (reads_left, nbadreads)

def pileup_to_msa(pileup):
    lines = 0
    columns = []
    noise = list()
    td_reads_set = set() #Store here all reads analysed
    nd_reads_set = set() #Store here all reads analysed
    td_msa_dict = dict() #Store here the allele of each read in each position
    nd_msa_dict = dict() #Store here the allele of each read in each position

    for line in pileup:
        column = line.strip().split()
        chrom, pos, ref = column[0], column[1], column[2]
        coverage_td, variants_td, reads_td = int(column[3]), column[4].upper(), column[6].split(',')
        coverage_nd, variants_nd, reads_nd = int(column[7]), column[8].upper(), column[10].split(',')

        ## Correct the variants lists and get the indels for the variants list ##
        corrected_variants_td, indel_list_td = correct_variants_list(variants_td)
        corrected_variants_nd, indel_list_nd = correct_variants_list(variants_nd)
        noise.append(len([i for i, x in enumerate(corrected_variants_td) if x not in [",", "."]]))

        ##Create tumor msa##
        done_reads = set()
        for index in range(0, len(reads_td)):
            if reads_td[index] not in td_msa_dict:
                td_msa_dict[reads_td[index]] = ["*"] * lines # Create a list for each read with as many asteriscs as positions in which it didn't appear
                td_msa_dict[reads_td[index]].append(corrected_variants_td[index]) #Append the variant
                done_reads.add(reads_td[index])
            elif reads_td[index] not in done_reads: #Make sure that we haven't read that read. Sometimes both mates overlap.
                    td_msa_dict[reads_td[index]].append(corrected_variants_td[index]) #Append the variant
                    done_reads.add(reads_td[index])

        ##Create control msa##
        done_reads = set()
        for index in range(0, len(reads_nd)):
            if reads_nd[index] not in nd_msa_dict:
                nd_msa_dict[reads_nd[index]] = ["*"] * lines # Create a list for each read with as many asteriscs as positions in which it didn't appear
                nd_msa_dict[reads_nd[index]].append(corrected_variants_nd[index]) #Append the variant
                done_reads.add(reads_nd[index])
            elif reads_nd[index] not in done_reads: #Make sure that we haven't read that read. Sometimes both mates overlap.
                nd_msa_dict[reads_nd[index]].append(corrected_variants_nd[index]) #Append the variant
                done_reads.add(reads_nd[index])

        td_reads_set = td_reads_set | set(reads_td) #Append all reads to the set
        nd_reads_set = nd_reads_set | set(reads_nd) #Append all reads to the set

        td_to_complete = td_reads_set - set(reads_td)
        for each in td_to_complete:
            td_msa_dict[each].append("*")

        nd_to_complete = nd_reads_set - set(reads_nd)
        for each in nd_to_complete:
            nd_msa_dict[each].append("*")

        lines +=1
        columns.append(pos)
    ## Do the msa ##
    td_msa = pd.DataFrame(td_msa_dict).T
    td_msa.columns = columns

    nd_msa = pd.DataFrame(nd_msa_dict).T
    nd_msa.columns = columns

    if "*" in td_msa.index: #Remove "*" as a row. It appears when in some positions there aren't reads and the pileup displays an asterisc.
        td_msa.drop("*", axis=0, inplace=True)

    if "*" in nd_msa.index: #Remove "*" as a row. It appears when in some positions there aren't reads and the pileup displays an asterisc.
        nd_msa.drop("*", axis=0, inplace=True)

    mean_noise = int(statistics.mean(noise))
    stdev_noise = int(statistics.stdev(noise))

    return(td_msa, nd_msa, mean_noise, stdev_noise)

def extract_context(msa, mut_pos, mut_base):
    context = set()
    context_frq = dict()
    bad_positions = 0
    seq_errors = 0
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
                seq_errors += 1
                continue

            if context_discrepancies == 0:
                context.add(pos+"_"+variant)
                context_frq[pos+"_"+variant] = variant_actual_frq
            elif context_discrepancies == 1: #One discrepant read may be due just a sequencing error
                context.add(pos+"_"+variant)
                context_frq[pos+"_"+variant] = variant_actual_frq
                seq_errors += 1
            else:
                bad_positions += 1
    return(context, bad_positions, seq_errors, context_frq)

def find_non_mutant_context(msa, context, mut_pos, mut_base):
    splitted_context = dict() #Store the context elements as pos:nucleotide
    for each in context:
        pos, nucleotide = each.split("_")
        splitted_context[pos] = nucleotide
    #Filter the matrix and keep non mutated reads that map at the mut position
    filtered_msa = msa[(msa[str(mut_pos)] != mut_base) & (msa[str(mut_pos)] != "*")] #Don't count the mutant reads. We want to know if context exist apart from the mutation. Also remove those that don't reach the mutation
    for each in splitted_context.keys():
        filtered_msa = filtered_msa[(filtered_msa[str(each)] == splitted_context[each]) | (filtered_msa[str(each)] == '*')]
    #Remove reads that contain too many asteriscs
    ctxt_msa = filtered_msa[splitted_context.keys()].replace("*", np.nan)
    min_length_coincidence = int(len(context) * 0.6) #Only allow two asteriscs per five contenxt elements
    ctxt_msa = ctxt_msa.dropna(thresh = min_length_coincidence)
    #Count mut reads with the context
    filtered_msa = msa[(msa[str(mut_pos)] == mut_base)] #Don't count the mutant reads. We want to know if context exist apart from the mutation. Also remove those that don't reach the mutation
    for each in splitted_context.keys():
        filtered_msa = filtered_msa[(filtered_msa[str(each)] == splitted_context[each]) | (filtered_msa[str(each)] == '*')]
    #Remove reads that contain too many asteriscs
    mut_msa = filtered_msa[splitted_context.keys()].replace("*", np.nan)
    mut_msa = mut_msa.dropna(thresh = min_length_coincidence)
    return (len(ctxt_msa), len(mut_msa))

def analyse_context(chrom, mut_pos, mut_base):
    ## Extract the pileup ##
    pileup_start = mut_pos - args.read_length if mut_pos >= args.read_length else 1 #Start can't be negative
    pileup_end = mut_pos + args.read_length
    pileup_command = ["samtools", "mpileup", "--output-QNAME", "-q", "30", "-Q", str(args.base_quality), "-R", "-f",  args.reference_fasta, args.tumor_bam, args.control_bam, "-r", chrom+":"+str(pileup_start)+"-"+str(pileup_end)]
    pileup = check_output(pileup_command).decode().splitlines() #run samtools mpileup

    ## Convert the pileup into a msa ##
    td_msa, nd_msa, mean_noise, stdev_noise = pileup_to_msa(pileup)

    ## Analyse the tumor msa to extract the context ##
    context, bad_positions, seq_errors, context_frq = extract_context(td_msa, mut_pos, mut_base)

    if bad_positions > 0.1 * len(context): #If greater, we'll discard de variant.
        tumor_non_mut_context_length = 0 #Set values for tumor_non_mut_context_length and control_non_mut_context_length so we can end the function without crashing.
        control_non_mut_context_length = 0
    else:    # Find the context in the control sample and in unmutated tumoral reads
        tumor_non_mut_context_length, tumor_mut_reads = find_non_mutant_context(td_msa, context, mut_pos, mut_base)
        control_non_mut_context_length, control_mut_reads = find_non_mutant_context(nd_msa, context, mut_pos, mut_base)

    return(tumor_mut_reads, tumor_non_mut_context_length, control_non_mut_context_length, control_mut_reads, len(context), seq_errors, mean_noise, stdev_noise)

def main_function(line):
    column = line.strip().split()
    chrom, pos, ref = column[0], column[1], column[2]
    coverage_td, variants_td, qual_td, reads_td = int(column[3]), column[4].upper(), column[5], column[6].split(',')
    coverage_nd, variants_nd = int(column[7]), column[8].upper()

    ##Skip flanking regions
    pos = int(pos)
    st, end = chrom.split(":")[1].split("-")
    if pos < 95 or pos > int(end) - int(st) - 95: # We allow 5 bp because it's still close and maybe we find some splicing mutations
        return
    ## Filter by minimun control coverage
    expected_nd_cov = coverage_td*args.control_coverage/100
    if coverage_nd < expected_nd_cov or coverage_td < args.tumor_coverage * 1.5: #If the genome is at 30x, at least we require 45x in each position as it should be repetitive
        if args.full:
            string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+"*"+"\t"+"*"
            if coverage_td < args.tumor_coverage * 1.5:
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

    real_alts_td, bad_alts = most_common_variant(variants_list_td_set, variants_list_td, variants_list_nd, args.tumor_threshold, args.control_max, coverage_td, coverage_nd)

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
        if len(element[1]) > 1 or len(element[0]) > 1: #Don't consider indels
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "indel")
            else:
                continue
        else:
            pass

        ##Launch polyfilter
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
        tumor_mut_reads, tumor_non_mut_context_length, control_non_mut_context_length, control_mut_reads, context_length, seq_errors, mean_noise, stdev_noise = analyse_context(chrom, pos, mut_base)

        if context_length != 0 and (tumor_non_mut_context_length < 1 or control_non_mut_context_length < 1) and (tumor_non_mut_context_length < 3 and control_non_mut_context_length < 3):
            if args.full: #Remove the mutation if the context only exists in the mutant reads
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "context_issues")
                continue
            else:
                continue
        else:
            pass

        if control_mut_reads > args.control_max:
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "germline_change")
                continue
            else:
                continue

        if tumor_mut_reads < args.tumor_threshold:
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "not_enough_reads")
                continue
            else:
                continue

        control_mutcov = variants_list_nd.count(mut_base)
        characteristics = [len(reads_left), nbadreads, seq_errors, coverage_td, control_mutcov, coverage_nd, mean_noise, stdev_noise]
        characteristics = list(map(str, characteristics))
        print(chrom, pos, args.name, element[0], element[1], ','.join(characteristics), ','.join(reads_left),"\n", sep = '\t', end = '')

if __name__ == '__main__':
    args = parse_args()

    arg_list = list()
    for arg in sorted(vars(args)):
        arg_list.append('--'+arg)
        arg_list.append(str(getattr(args, arg)))
    ## VCF header
    print('##Command = python3 %s %s' % (argv[0], ' '.join(arg_list)))

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
