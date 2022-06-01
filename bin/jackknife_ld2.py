#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description='Calculate the jacknife estimates of LD')
parser.add_argument('bfile', help='Input bfile prefix')
parser.add_argument('keep', help='Keep samples')
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
  dtype='string', names=['FID', 'IID'], usecols=[0,1]
)

# join on iids
iid = pd.DataFrame(snp_on_disk.iid, columns=['FID', 'IID'])
iid = iid.set_index(['FID', 'IID'])
keep = keep.join(iid, ['FID', 'IID'], 'inner')

# get the indices
iid_index = snp_on_disk.iid_to_index(keep.values)

# read the data
snpdata = snp_on_disk[iid_index,:].read().standardize()

# jacknife function
def cor_jackknife(A, B):
    A2 = A**2
    B2 = B**2
    AB = A * B[:,None]

    Asum = A.sum(axis=0)
    Bsum = B.sum()

    A2sum = A2.sum(axis=0)
    B2sum = B2.sum()

    ABsum = AB.sum(axis=0)

    Aisum   = Asum[None, :] - A
    Bisum   = Bsum - B

    A2isum  = A2sum[None, :] - A2
    B2isum  = B2sum - B2

    ABisum = ABsum[None, :] - AB

    N = len(B)

    r = (N*ABsum - Asum * Bsum) / np.sqrt((N*A2sum - Asum**2) * (N*B2sum - Bsum**2))
    ris = ((N-1)*ABisum - Aisum*Bisum[:,None]) / np.sqrt(((N-1)*A2isum - Aisum**2) * ((N-1)*B2isum - Bisum**2)[:,None])
    rimean = ris.mean(axis=0)
    ri2mean = (ris**2).mean(axis=0)

    return(N*r - (N-1)*rimean, N*(r**2) - (N-1)*ri2mean)

# set up r and r2
M = snpdata.sid_count
N = int(M*(M+1)/2)

rs  = np.ones(N)
r2s = np.ones(N)

k = 1
for i in range(1,M):
    r, r2 = cor_jackknife(snpdata.val[:,:i], snpdata.	val[:,i])

    rs[k:(k+i)] = r
    r2s[k:(k+i)] = r2

    k += i + 1

# write it out
rs.tofile(args.outr)
r2s.tofile(args.outr2)
