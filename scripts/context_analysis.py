#!/usr/bin/env python3
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
            filtered_msa = filtered_msa[(filtered_msa[str(each)] == splitted_context[each]) | (filtered_msa[str(each)] == '*')]
        
        #Remove reads that contain too many asteriscs
        mut_msa = filtered_msa[splitted_context.keys()].replace("*", np.nan)
        mut_msa = mut_msa.dropna(thresh = min_length_coincidence)
        mut_names = mut_msa.index.values.tolist()  
        full_mut_msa = msa.loc[mut_names]

    #Check if there's strand bias
    posteriors_strand = bayes_strand.strand_bias(full_mut_msa, full_ctxt_msa)
    
    return (ctxt_length, len(full_mut_msa), posteriors_strand)

def analyse_context(chrom, mut_pos, mut_base):
    ## Extract the pileup ##    
    tumor_mq = get_pileup(tumor_bam, chrom, mut_pos, fasta)
    control_mq = get_pileup(control_bam, chrom, mut_pos, fasta)
    
    ## Convert the pileup into a msa ##
    td_msa, mean_noise, stdev_noise = pileup_to_msa(tumor_mq)
    nd_msa = pileup_to_msa(control_mq)[0]

    ## Analyse the tumor msa to extract the context ##
    context, bad_positions, seq_errors, context_frq = extract_context(td_msa, mut_pos, mut_base)
    
    if bad_positions > 0.1 * len(context): #If greater, we'll discard de variant. Make all returning variables null
        tumor_non_mut_context_length, tumor_mut_reads, control_mut_reads, control_non_mut_context_length, posteriors_strand_tumor, posteriors_strand_control = 0, 0, 0, 0, (0,0), (0,0)
    else:    # Find the context in the control sample and in unmutated tumoral reads
        tumor_non_mut_context_length, tumor_mut_reads, posteriors_strand_tumor = find_non_mutant_context(td_msa, context, mut_pos, mut_base)
        control_non_mut_context_length, control_mut_reads, posteriors_strand_control = find_non_mutant_context(nd_msa, context, mut_pos, mut_base)

    return(tumor_mut_reads, tumor_non_mut_context_length, control_non_mut_context_length, control_mut_reads, len(context), seq_errors, mean_noise, stdev_noise, posteriors_strand_tumor)
