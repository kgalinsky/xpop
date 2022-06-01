#!/usr/bin/env python

import argparse

import numpy as np
import numpy.random as nr
import scipy.stats as ss

parser = argparse.ArgumentParser(description='Generate simulation.')
parser.add_argument('ldfiles', type=argparse.FileType('rb', 0), nargs=2,
                   help='LD files')
parser.add_argument('-l, --lambda', metavar='l', type=float, default=0
                   help='lambda for ridge regression')

args = parser.parse_args()

print args