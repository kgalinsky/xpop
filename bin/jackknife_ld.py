#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description='Calculate the jacknife estimates of LD')
parser.add_argument('bfile', help='Input bfile prefix')
parser.add_argument('keep', help='Keep SNPs')
parser.add_argument('outr', help='Output r')
parser.add_argument('outr2', help='Output r^2')
args = parser.parse_args()

import pandas as pd
import numpy as np
from pysnptools.snpreader import Bed

# open everything up
snp_on_disk = Bed(args.bfile, count_A1=True)

keep = pd.read_table(
  args.keep, delim_whitespace=True,
  dtype='string', names=['FID', 'IID', 'CLST'], index_col='CLST'
)

# join on iids
iid = pd.DataFrame(snp_on_disk.iid, columns=['FID', 'IID'])
iid = iid.set_index(['FID', 'IID'])
keep = keep.join(iid, ['FID', 'IID'], 'inner')

# get the indices
iid_index = snp_on_disk.iid_to_index(keep.values)

# read the data
snpdata = snp_on_disk[iid_index,:].read()

# jacknife function
def cor_jacknife(X):
    X = X[~np.isnan(X).any(axis=1)]

    X2  = X**2
    X12 = X.prod(axis=1)

    Xsum   = X.sum(axis=0)
    X2sum  = X2.sum(axis=0)
    X12sum = X12.sum()
    
    Xisum   = Xsum[None, :] - X
    X2isum  = X2sum[None, :] - X2
    X12isum = X12sum - X12
    
    N = len(X12)
    
    r = (N*X12sum - Xsum.prod()) / np.sqrt((N*X2sum - Xsum**2).prod())
    ris = ((N-1)*X12isum - Xisum.prod(axis=1)) / np.sqrt(((N-1)*X2isum - Xisum**2).prod(axis=1))
    rimean = ris.mean()
    ri2mean = (ris**2).mean()
    
    return(N*r - (N-1)*rimean, N*(r**2) - (N-1)*ri2mean)

# set up r and r2
M = snpdata.sid_count
N = M*(M+1)/2
rs  = np.empty(N)
r2s = np.empty(N)

k = 0
for i in range(M):
    for j in range(i):
        r, r2 = cor_jacknife(snpdata.val[:,(i,j)])
        rs[k] = r
        r2s[k] = r2
        k += 1
    rs[k] = 1
    r2s[k] = 1
    k = k + 1

# write it out
rs.tofile(args.outr)
r2s.tofile(args.outr2)
