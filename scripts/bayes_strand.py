#!/usr/bin/env python3
from scipy.stats import beta

def contingency_bayes(x1, y1, x2, y2, size = 10000, tail = 'both', seed = 1234): 
    tail = tail.lower()
    theta1 = beta.rvs(x1+1, y1+1, size = size, random_state = seed) 
    theta2 = beta.rvs(x2+1, y2+1, size = size, random_state = seed)
    if tail == 'lower':
        result = sum(theta1<theta2)/size
    elif tail == 'greater':
        result = sum(theta1>theta2)/size
    elif tail == 'both':
        result = [sum(theta1>theta2)/size, sum(theta1<theta2)/size]
    else:
        raise ValueError("'tail' argument accepted values are: 'lower', 'greater' or 'both'")
    return result

def msa_strand_count(msa):
    fwd = 0
    rev = 0
    if len(msa) > 0:
        for index, rows in msa.iterrows(): 
            rows = [str(w).replace(" ", "") for w in rows]
            f = rows.count(",")
            r = rows.count(".")
            #assert(any(i == 0 for i in [f, r])), "%s read is both strand and reverse. %s" % (index, rows)
            if f > r:
                fwd += 1
            elif r > f: # We won't find r == f, but if it happened, we wouldn't count that read
                rev += 1
    else:
        pass
    return fwd, rev

def strand_bias(msa1, msa2):
    fwd1, rev1 = msa_strand_count(msa1)
    fwd2, rev2 = msa_strand_count(msa2)
    prev, pfwd = contingency_bayes(rev1, fwd1, rev2, fwd2, 10000, 'both')
    return prev, pfwd
