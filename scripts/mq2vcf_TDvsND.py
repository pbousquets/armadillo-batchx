#!/usr/bin/env python3
from pickle import TRUE
from sys import argv, stdin
from re import findall, match
from multiprocessing import Pool
from subprocess import check_output
from fit import fitArmNet
import argparse
import pysam 
import statistics
import pandas as pd
import numpy as np
import bayes_strand 


def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = argparse.ArgumentParser(description = 'Converts an MQ file into a VCF filtering by coverage')
    parser.add_argument(
    '-cc', '--control_coverage', type = int, required = False, default = 30, metavar = 'INT',
    help = 'Control genome coverage (default: %(default)s)')
    parser.add_argument(
    '-cb', '--control_bam', type = str, required = True, metavar = 'FILE',
    help = 'Tumor miniBAM to analyse')
    parser.add_argument(
    '-tc', '--tumor_coverage', type = int, required = False, default = 30, metavar = 'INT',
    help = 'Tumor genome coverage')
    parser.add_argument(
    '-f', '--full', action = 'store_true', default = True,
    help = 'Print log with discarded variants and add a mutant readname column.')
    parser.add_argument(
    '-i', '--input', type = str, required = False, default = '-', metavar = 'FILE',
    help = 'MQ file to convert')
    parser.add_argument(
    '-gc', '--GCcutoff', type = int, required = False, default = 80, metavar = 'INT',
    help = 'Max GC% per read. (default: %(default)s)')
    parser.add_argument(
    '-m', '--model', type = str, required = False, default = 'armadillo/lib/armNet_epoch80.pth',
    help = 'pretrained CNN model to fit the mutations. (default: %(default)s)')
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
    '-Q', '--base_quality', type = int, required = False, default = 30,
    help = 'Filter by base quality value (default: %(default)s)')
    parser.add_argument(
    '-q', '--mapq', type = int, required = False, default = 30,
    help = 'Filter by mapping quality (default: %(default)s)')    
    parser.add_argument(
    '-r', '--reference_fasta', type = str, required = True,
    help = 'Reference fasta file')
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

def most_common_variant(single_variants_list, full_variants_list_TD, full_variants_list_ND, threshold_TD, coverage_td, coverage_nd):
    alts = list()
    germline_list = list()
    overlimit = list()
    for variant in single_variants_list:
        count_TD = full_variants_list_TD.count(variant)
        count_ND = full_variants_list_ND.count(variant)
        if 1.5 * args.tumor_coverage > count_TD >= threshold_TD and count_ND <= args.control_max: 
            pp = bayes_strand.contingency_bayes(count_TD, coverage_td, count_ND, coverage_nd, 10000, "greater")
            if pp < 0.9: 
                germline_list.append(variant)  #Store a list of those that were discarded due to high freq in control (likely SNPs and quimeric variants due to the mapping step)
            else:
                alts.append(variant)
        elif 1.5 * args.tumor_coverage > count_TD >= threshold_TD and count_ND > args.control_max:
            germline_list.append(variant)
        elif 1.5 * args.tumor_coverage < count_TD >= threshold_TD and count_ND <= args.control_max:
            overlimit.append(variant)
        else:
            pass
    return(alts, germline_list, overlimit)

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

def view2fasta(bamdir, coord, reads): # Samtools: get reads' sequences
    samtools_command = ['samtools', 'view', bamdir, coord ]
    fastareads = ''
    result = check_output(samtools_command).decode().splitlines()
    for readres in result:
        column = readres.split("\t")
        if column[0] in reads: #If the read is in the list of reads to analyse and isn't bad, analyse it
            GC = GC_filt(column[9]) #GC filter. If excesive GC content, remove the read
            if GC <= args.GCcutoff:
                fastareads = fastareads+">"+column[0]+"\n"+column[9]+"\n" #Store as fasta.
    return fastareads

def blat_filter(blat_result): #Filter the blat result to remove the reads with 100% identity
    badreads, allreads, reads_left = list(), list(), list()
    for line in blat_result:
        column = line.strip().split("\t")
        try:
            ID = column[9]
            if ID in badreads:
                continue
            if ID not in allreads:
                allreads.append(ID) #store all IDs in a list
            if column[0] == column[10] and all(int(i) == 0 for i in column[1:8]) and int(column[11]) == 0 and column[12] == column[10] and ID not in badreads: #if perfect match, the change in the read is not a mutation..
                badreads.append(ID) #.. so append it to the bad reads list
                continue
        except IndexError:
            pass
    reads_left = list(set(allreads) - set(badreads))
    return(reads_left)

def get_pileup(bam, chr, mutpos, reference_genome, mapq, bq):
    mutpos = int(mutpos)
    pileup = list()
    for column in bam.pileup(chr, mutpos, mutpos+1, stepper = "samtools", fastafile = reference_genome, min_base_quality = bq, min_mapping_quality = mapq):
        pos = column.reference_pos + 1
        alts = column.get_query_sequences(mark_matches=True, mark_ends=False, add_indels=True)
        reads = column.get_query_names()
        pileup.append([pos, alts, reads])
    return pileup

def pileup_to_msa(pileup):
    lines = 0
    columns = []
    noise = list()
    reads_set = set() #Store here all reads analysed
    msa_dict = dict() #Store here the allele of each read in each position

    for element in pileup:
        pos, variants, reads = element[0], element[1], element[2]

        ## Correct the variants lists and get the indels for the variants list ##
        corrected_variants, _ = correct_variants_list(variants)
        noise.append(len([i for i, x in enumerate(corrected_variants) if x not in [",", "."]]))

        ##Create msa##
        done_reads =set()
        for index in range(0, len(reads)):
            if reads[index] not in msa_dict:
                msa_dict[reads[index]] = ["*"] * lines # Create a list for each read with as many asteriscs as positions in which it didn't appear
                msa_dict[reads[index]].append(corrected_variants[index]) #Append the variant
                done_reads.add(reads[index])
            elif reads[index] not in done_reads: #Make sure that we haven't read that read. Sometimes both mates overlap.
                    msa_dict[reads[index]].append(corrected_variants[index]) #Append the variant
                    done_reads.add(reads[index])

        reads_set = reads_set | set(reads) #Append all reads to the set

        to_complete = reads_set - set(reads)
        for each in to_complete:
            msa_dict[each].append("*")

        lines +=1
        columns.append(str(pos))
        
    ## Get the msa ##
    msa = pd.DataFrame(msa_dict).T
    msa.columns = columns

    if "*" in msa.index: #Remove "*" as a row. It appears when in some positions there aren't reads and the pileup displays an asterisc.
        msa.drop("*", axis=0, inplace=True)
    msa = msa.applymap(lambda s:s.upper() if type(s) == str else s)

    mean_noise = int(statistics.mean(noise))
    stdev_noise = int(statistics.stdev(noise))

    return(msa, mean_noise, stdev_noise)
    
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
            variant_frq = only_mut_msa[pos].tolist().count(variant) 
            context_discrepancies = len(only_mut_msa[pos].tolist()) - only_mut_msa[pos].tolist().count("*") - variant_frq #"*" means the read doesn't map there, so we cannot count that one.
            if variant_frq == 1 and context_discrepancies > 0: #if the variant appears just once, it may be just a sequencing error
                seq_errors += 1
                continue

            if context_discrepancies == 0:
                context.add(pos+"_"+variant)
                context_frq[pos+"_"+variant] = variant_frq
            elif context_discrepancies == 1: #One discrepant read may be due just a sequencing error
                context.add(pos+"_"+variant)
                context_frq[pos+"_"+variant] = variant_frq
                seq_errors += 1
            else:
                bad_positions += 1
    return(context, bad_positions, seq_errors, context_frq)

def remove_other_contexts(df):
    df = df.applymap(lambda s:s.upper() if type(s) == str else s) #Convert bases to uppercase
    df_counts = df.apply(pd.Series.value_counts) #Summarize
    counts = df_counts[df_counts.index.isin(['A', 'T', 'C', 'G'])] #Filter by ATCG nucleotides
    rm_cols = counts[counts.columns[counts.max() > 5]].dropna(how='all') #Keep the consistent ones that could define different haplotypes
    
    #Get a dictionary with pos - nucleotides to filter de df
    dic = {}
    for col_idx, col in rm_cols.iteritems():
        na = col[~np.isnan(col)]
        for each in na.index:
            if na[each] > 5 and col_idx not in dic.keys():
                dic[col_idx] = [each]
            elif col_idx in dic.keys():
                dic[col_idx].append(each)
                
    #Filter the df          
    for key, value in dic.items():
        for each in value:
            df = df[df[key] != each]
    return df

def find_non_mutant_context(msa, context, mut_pos, mut_base):
    if len(context) == 0:
        full_ctxt_msa = msa[(msa[str(mut_pos)] != mut_base)]
        full_ctxt_msa = remove_other_contexts(full_ctxt_msa)
        ctxt_length = len(full_ctxt_msa)
        full_mut_msa = msa[(msa[str(mut_pos)] == mut_base)]
    else:
        splitted_context = dict()
        for each in context:
            pos, nucleotide = each.split("_")
            if str(each.split("_")[0]) in list(msa.columns.values):
                splitted_context[pos] = nucleotide

        #Filter the matrix and keep non mutated reads that map at the mut position
        filtered_msa = msa[(msa[str(mut_pos)] != mut_base) & (msa[str(mut_pos)] != "*")] #Don't count the mutant reads. We want to know if context exist apart from the mutation. Also remove those that don't reach the mutation
        for each in splitted_context.keys():
            filtered_msa = filtered_msa[(filtered_msa[str(each)] == splitted_context[each]) | (filtered_msa[str(each)] == '*')]
                
        #Remove reads that contain too many asteriscs
        ctxt_msa = filtered_msa[splitted_context.keys()].replace("*", np.nan)
        min_length_coincidence = int(len(context) * 0.75) #Only allow one asterisc per four context elements
        ctxt_msa = ctxt_msa.dropna(thresh = min_length_coincidence)
        unmut_names = ctxt_msa.index.values.tolist()  
        full_ctxt_msa = msa.loc[unmut_names]
        ctxt_length = len(full_ctxt_msa)
        column_nas = [ctxt_length - ctxt_msa[col].count() for col in ctxt_msa.columns]
        if any(excesive_na > ctxt_length * 0.5 for excesive_na in column_nas) and len(ctxt_msa.dropna(thresh = min_length_coincidence)) < 3:
            full_ctxt_msa = [] # If the context isn't perfect in at least 3 reads and any column isn't well represented (more than the half reads contain NAs), remove the matrix to make sure that the mutation is discarded afterwards
        else:
            pass
        
        #Count mut reads with the context
        filtered_msa = msa[(msa[str(mut_pos)] == mut_base)] 
        for each in splitted_context.keys():
            if str(each) in filtered_msa.index:
                filtered_msa = filtered_msa[(filtered_msa[str(each)] == splitted_context[each]) | (filtered_msa[str(each)] == '*')]

        #Remove reads that contain too many asteriscs
        mut_msa = filtered_msa
        mut_names = mut_msa.index.values.tolist()  
        full_mut_msa = msa.loc[mut_names]

    #Check if there's strand bias
    posteriors_strand = bayes_strand.strand_bias(full_mut_msa, full_ctxt_msa)
    return (ctxt_length, full_mut_msa.index.tolist(), posteriors_strand)

def analyse_context(chrom, mut_pos, mut_base, fasta, tumor_bam, control_bam):
    ## Extract the pileup ##    
    tumor_mq = get_pileup(tumor_bam, chrom, mut_pos, fasta, args.mapq, args.base_quality)
    control_mq = get_pileup(control_bam, chrom, mut_pos, fasta, args.mapq -10, args.base_quality -10)
    
    ## Convert the pileup into a msa ##
    td_msa, mean_noise, stdev_noise = pileup_to_msa(tumor_mq)
    nd_msa = pileup_to_msa(control_mq)[0]

    ## Analyse the tumor msa to extract the context ##
    context, bad_positions, seq_errors, context_frq = extract_context(td_msa, mut_pos, mut_base)
    
    if bad_positions > 0.1 * len(context): #If greater, we'll discard de variant. Make all returning variables null
        tumor_non_mut_context_length, tumor_mut_reads, control_mut_reads, control_non_mut_context_length, posteriors_strand_tumor, posteriors_strand_control = 0, [], [], 0, (0,0), (0,0)
    else:    # Find the context in the control sample and in unmutated tumoral reads
        tumor_non_mut_context_length, tumor_mut_reads, posteriors_strand_tumor = find_non_mutant_context(td_msa, context, mut_pos, mut_base)
        control_non_mut_context_length, control_mut_reads, posteriors_strand_control = find_non_mutant_context(nd_msa, context, mut_pos, mut_base)

    return(tumor_mut_reads, tumor_non_mut_context_length, control_non_mut_context_length, control_mut_reads, len(context), seq_errors, mean_noise, stdev_noise, posteriors_strand_tumor)

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
    
    ## Filter by minimum control coverage
    if coverage_nd < args.control_coverage * 0.8 or coverage_td < args.tumor_coverage * 1.5: #If the genome is at 30x, at least we require 45x in each position as it should be repetitive
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

    real_alts_td, germline_list, overlimit_list = most_common_variant(variants_list_td_set, variants_list_td, variants_list_nd, args.tumor_threshold, coverage_td, coverage_nd)
       
    ## Annotate alts with frequency over threshold in control sample.
    if args.full:
        for alt in germline_list:
            string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref+"\t"+alt
            print_log(string, "germline_change")
        for alt in overlimit_list:
            string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref+"\t"+alt
            print_log(string, "VAF_over_coverage")
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

        ##Launch blat
        coord = chrom+":"+str(pos)+"-"+str(pos)
        reads_fasta = view2fasta(args.tumor_bam, coord, real_reads_td)

        if len(reads_fasta) == 0: #If no reads left, stop
            continue
        else:
            pass

        reads_left = reads_fasta
        if len(reads_left) < args.tumor_threshold: #If there are enough reads, go on
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+element[0]+"\t"+element[1]
                print_log(string, "blat_perfect_match")
                continue
            else:
                continue
        else:
            pass

        ## Phase the variants to check if all the reads containing the mutation come from the same region
        fasta = pysam.FastaFile(args.reference_fasta)
        tumor_bam = pysam.AlignmentFile(args.tumor_bam, "rb" )
        control_bam = pysam.AlignmentFile(args.control_bam, "rb" )
        ref_base, mut_base = element[0], element[1]
        tumor_mut_reads, tumor_non_mut_context_length, control_non_mut_context_length, control_mut_reads, context_length, seq_errors, mean_noise, stdev_noise, posteriors_strand_tumor = analyse_context(chrom, pos, mut_base, fasta, tumor_bam, control_bam)
        tumor_mut_reads_len, control_mut_reads_len = len(tumor_mut_reads), len(control_mut_reads)
        pp = bayes_strand.contingency_bayes(tumor_mut_reads_len, tumor_non_mut_context_length, control_mut_reads_len, control_non_mut_context_length, 10000, "greater")
        rbias, fbias = posteriors_strand_tumor

        ## Last filters
        if tumor_mut_reads_len < args.tumor_threshold:
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref_base+"\t"+mut_base
                print_log(string, "not_enough_reads")
                continue
            else:
                continue
        
        if control_mut_reads_len > args.control_max:
            if args.full:
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref_base+"\t"+mut_base
                print_log(string, "germline_change")
                continue
            else:
                continue

        if context_length != 0 and (tumor_non_mut_context_length + tumor_mut_reads_len < args.tumor_coverage * 0.4 or control_non_mut_context_length + control_mut_reads_len < args.control_coverage * 0.75):
            if args.full: #Remove the mutation if the context only exists in the mutant reads
                string = chrom+"\t"+str(pos)+"\t"+args.name+"\t"+ref_base+"\t"+mut_base
                print_log(string, "low_context_coverage")
                continue
            else:
                continue
        else:
            pass
        
        ## Fit to CNN
        keep, qual = fitArmNet(chrom, int(pos), mut_base, args.tumor_bam, args.control_bam, args.reference_fasta, args.model)
        
        if pp < 0.95:
            FILTER = "GERM"
        else:
            FILTER = "PASS"

        FILTER = FILTER if keep else "POTENTIAL_ARTIFACT" # Check if the model passed the mutation 

        unf_TD_count = [item.upper() for item in variants_list_td].count(mut_base)
        unf_ND_count = [item.upper() for item in variants_list_nd].count(mut_base)
        
        INFO = ["PP="+str(pp), "FB="+str(fbias), "RB="+str(rbias), "MN="+str(mean_noise), "SN="+str(stdev_noise)]
        FORMAT = "AD:uAD:CD:DP"
        TUMOR = [tumor_mut_reads_len, unf_TD_count, tumor_non_mut_context_length, coverage_td]
        CONTROL = [control_mut_reads_len, unf_ND_count, control_non_mut_context_length, coverage_nd]
        TUMOR, CONTROL = list(map(str, TUMOR)), list(map(str, CONTROL))
        
        print_list = [str(chrom), str(pos), args.name, element[0], mut_base, str(qual), FILTER, ";".join(INFO), FORMAT, ":".join(TUMOR), ":".join(CONTROL), ",".join(tumor_mut_reads)]

        return print_list

if __name__ == '__main__':
    args = parse_args()
    f = open(f'{args.name}_candidates.vcf', "w+")
    ## VCF header
    arg_list = list()
    for arg in sorted(vars(args)):
        arg_list.append('--'+arg)
        arg_list.append(str(getattr(args, arg)))

    print('##fileformat=VCFv4.2',
    '##FILTER=<ID=PASS,Description="All filters passed">',
    '##FILTER=<ID=GERM,Description="Somatic posterior probability < 0.95, potential germline event">',
    '##FILTER=<ID=POTENTIAL_ARTIFACT,Description="Classified as artifact by the model">',
    '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Count of mutant reads">',
    '##FORMAT=<ID=uAD,Number=1,Type=Integer,Description="Count of unfiltered mutant reads">',
    '##FORMAT=<ID=CD,Number=1,Type=Integer,Description="Count of non-mutant context reads">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Count of total reads">',
    '##INFO=<ID=PP,Number=1,Type=Float,Description="Somatic mutation posterior probability">',
    '##INFO=<ID=FB,Number=1,Type=Float,Description="Forward strand bias posterior probability">',
    '##INFO=<ID=RB,Number=1,Type=Float,Description="Reverse strand bias posterior probability">',
    '##INFO=<ID=MN,Number=1,Type=Float,Description="Mean alternative allele frequency per position in context">',
    '##INFO=<ID=NM,Number=1,Type=Float,Description="Standard deviation of alternative allele frequency per position in context">',
    '##Command = python3 %s %s' % (argv[0], ' '.join(arg_list)),
    "\t".join(['#CHROM','POS','ID','REF','ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.tumor_bam, args.control_bam]),  file = f, sep = "\n",)

    ## Using multiprocessing.Pool each task is run in a differente thread with their own variables, and the result is given only when all the task are completed.
    pool = Pool(processes = args.threads)
    ## Initialize the pool of threads
    if args.input == '-':
        results = pool.imap(func = main_function, iterable = stdin)
    else:
        pileup = open(args.input)
        results = pool.imap(func = main_function, iterable = pileup)
    pool.close()
    pool.join()
    results = filter(None, results)

    for mut in sorted(list(results)):
        print("\t".join(mut), file = f)
