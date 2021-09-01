import pandas as pd
import numpy as np
from torch import tensor, stack, device as torch_device
from torch.cuda import is_available
import pysam
from statistics import mean
from warnings import catch_warnings, simplefilter

class bam_features:
    def __init__(self, bamname, chrom, pos, fasta):
        self.bamname = bamname
        self.chrom = chrom
        self.pos = int(pos)
        self.bam = pysam.AlignmentFile(self.bamname, "rb")
        region = f"{chrom}:{pos-23}-{pos+26}"
        self.refseq = "".join(pysam.faidx(fasta, str(region)).strip().split("\n")[1:]) #Retrieve reference sequence
    
    def get_msa(self, reads, min_pos, max_pos):
        msa_dict = {"seqs": dict(), "quals": dict(), "poss": dict()}
        n=0
        for read in reads:
            n+=1
            id = f"{read.query_name}_{n}"
            query_seq = list(read.query_sequence) #includes softc
            query_q = list(read.query_qualities) #includes softc
            query_readpos = np.array(range(1, read.query_length +1))/read.query_length
            query_readpos = query_readpos.tolist()
            cigar = read.cigartuples
            left_softc = 0

            if len(cigar) > 1: #Check if there's softclipping on the left
                left = 0
                if cigar[0][0] == 4: #Left side
                    left_softc = cigar[0][1]

                for cig_type, length in cigar:
                    if cig_type == 1: #Insertion -> Convert the whole insertion to one single position to avoid misalignments
                        if left == 0: # Once I found a case with a starting insertion that crashed. Process slightly differently this insertions:
                            ins = "".join(query_seq[:length]).lower()
                            query_seq = [ins] + query_seq[left+length:]
                            query_q = [mean(query_q[:length])] + query_q[length:]
                            query_readpos = [mean(query_readpos[:length])] + query_readpos[length:]
                        else:
                            ins = "".join(query_seq[left-1:left+length]).lower()
                            query_seq = query_seq[:left-1] + [ins] + query_seq[left+length:]
                            query_q = query_q[:left-1] + [mean(query_q[left-1:left+length])] + query_q[left+length:]
                            query_readpos = query_readpos[:left-1] + [mean(query_readpos[left-1:left+length])] + query_readpos[left+length:]
                    elif cig_type == 2: #Deletion -> Fill in the gaps with "-"
                        query_seq = query_seq[:left] + ["-"] * length  + query_seq[left:]
                        query_q = query_q[:left] + [np.nan] * length + query_q[left:]
                        query_readpos = query_readpos[:left] + [np.nan] * length + query_readpos[left:]
                        left += length
                    else:
                        left += length
                    
            start = read.reference_start
            start_softc = start - left_softc
            positions = list(range(start_softc, start_softc + read.query_length + 1))

            if positions[-1] < max_pos: # If the read ends before the window ends, add NAs
                n = max_pos - positions[-1] + 1 
                positions += list(range(positions[-1] + 1, max_pos + 1))
                query_seq = query_seq + [np.nan] * n
                query_q = query_q + [np.nan] * n
                query_readpos = query_readpos + [np.nan] * n

            if start_softc > min_pos: # If the read starts after the window starts, add NAs
                n = start_softc - min_pos
                query_seq = [np.nan] * n + query_seq
                query_q = [np.nan] * n + query_q
                query_readpos = [np.nan] * n + query_readpos
                positions = list(range(min_pos, start_softc)) + positions

            st_idx = positions.index(min_pos)
            end_idx = positions.index(max_pos)+1

            query_seq = query_seq[st_idx:end_idx]
            query_q = query_q[st_idx:end_idx]
            query_readpos = query_readpos[st_idx:end_idx]

            msa_dict["seqs"][id] = list(query_seq)
            msa_dict["quals"][id] = query_q
            msa_dict["poss"][id] = query_readpos

        seqs = pd.DataFrame.from_dict(msa_dict["seqs"], orient='index', columns=list(range(min_pos+1, max_pos+2))) #Move 1 base due because we work with 1-based system
        quals = pd.DataFrame.from_dict(msa_dict["quals"], orient='index', columns=list(range(min_pos+1, max_pos+2)))
        poss = pd.DataFrame.from_dict(msa_dict["poss"], orient='index', columns=list(range(min_pos+1, max_pos+2)))
        return (seqs, quals, poss)

    def get_reads(self):
        min_pos = self.pos - 24
        max_pos = self.pos + 25
        self.reads = self.bam.fetch(self.chrom, min_pos, max_pos) #Different from get_reads function because one program is 0-based, the other 1-based
        return self.get_msa(self.reads, min_pos, max_pos)

def filter_msa_by_variants(msa, variants, filter_include = True, dropna = False):
    assert filter_include == True or filter_include == False, "filter_include expected to be a boolean"
    assert dropna == True or dropna == False, "dropna expected to be a boolean"

    for variant_pos, variant in variants.items():
        variant = [variant] if type(variant) is str else variant
        if filter_include:
            msa = msa[(msa[variant_pos].isin(variant)) | msa[variant_pos].isna()]
        else:
            msa = msa[~msa[variant_pos].isin(variant)]
    if dropna:
        msa = msa[msa[variant_pos].notna()]
    
    return msa

def extract_context(msa, ref):
    i = 0
    context = dict()
    for pos in msa.columns:
        refbase = ref[i]
        for variant in msa[pos].dropna().unique().tolist():
            if variant == refbase:
                continue
            variant_frq = msa[pos].tolist().count(variant) 
            if variant_frq > 1 and len(msa[pos].tolist()) - msa[pos].isna().sum():
                context[pos] = [variant] if pos not in context.keys() else context[pos] + [variant]
        i += 1            
        
    return(context)



def model_features(tumor_bam, control_bam, chrom, pos, alt, fasta):
    select = ["A", "T", "C", "G"]
    mut = {pos: alt}
    tumor = True
    res = list()
    device = torch_device("cpu")
    for bam in [tumor_bam, control_bam]:
        feat = bam_features(bam, chrom, pos, fasta)
        seqs, quals, poss = feat.get_reads() # MSA 
        mut_msa = filter_msa_by_variants(seqs, mut, filter_include=True, dropna = True) # Mut reads
        mreads = mut_msa.index.values.tolist()
        mut_quals = quals.loc[mreads] #quals of mut reads
        mut_poss = poss.loc[mreads] #positions of mut reads
        if tumor:
            ctxt = extract_context(mut_msa, feat.refseq)
            del ctxt[pos]

        wt_msa = filter_msa_by_variants(seqs, mut, filter_include=False, dropna=True)
        wt_msa = filter_msa_by_variants(wt_msa, ctxt, filter_include=True, dropna=False)

        creads = wt_msa[list(ctxt.keys())].index.values.tolist()
        wt_msa = wt_msa.loc[creads]
        wt_quals= quals.loc[creads] #quals of wt reads
        wt_poss = poss.loc[creads] #positions of wt reads

        ## Summarize ##
        with catch_warnings():
            simplefilter("ignore")
            if len(mreads) == 0:
                for each in select:
                    mut_msa.loc[each] = [0] * len(mut_msa.columns)
                    mut_quals.loc["quals"] = [0] * len(mut_msa.columns)
                    mut_poss.loc["poss"] = [0] * len(mut_msa.columns)
                mut_msa = mut_msa.infer_objects()
            else:
                mut_msa = mut_msa.fillna(value=-1).apply(pd.value_counts).loc[select].fillna(0)

            if len(creads) == 0:
                for each in select:
                    wt_msa.loc[each] = [0] * len(wt_msa.columns)
                    wt_quals.loc["quals"] = [0] * len(wt_msa.columns)
                    wt_poss.loc["poss"] = [0] * len(wt_msa.columns)
                wt_msa = wt_msa.infer_objects()
            else:
                wt_msa = wt_msa.fillna(value=-1).apply(pd.value_counts).loc[select].fillna(0)

            mut_quals = mut_quals.fillna(value=-1).median(axis = 0)
            mut_poss = mut_poss.fillna(value=-1).mean(axis = 0)
            wt_quals = wt_quals.fillna(value=-1).median(axis = 0)
            wt_poss = wt_poss.fillna(value=-1).mean(axis = 0)
        
        mut_msa.loc["q"] = mut_quals
        mut_msa.loc["poss"] = mut_poss
        wt_msa.loc["q"] = wt_quals
        wt_msa.loc["poss"] = wt_poss
        res += [tensor(np.array(mut_msa)).to(device), tensor(np.array(wt_msa)).to(device)]

        tumor = False

    ref = pd.DataFrame.from_dict({'ref': list(feat.refseq)}, orient='index', columns=wt_msa.columns.values.tolist()).apply(pd.value_counts).loc[select].fillna(0)
    ref.loc["mut"] = [0] * len(ref.columns)
    ref.at["mut", pos] = 1
    ref.loc["ctxt"] = [0] * len(ref.columns)
    ref.at["ctxt", ctxt.keys()] = 1    
    res += [tensor(np.array(ref)).to(device)]

    return stack(res)