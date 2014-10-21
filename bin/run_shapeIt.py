__author__ = 'Nicholas Harding'

# This script will submit a shapeIt run using the phasing package. Basically an
# alternative to the ipython notebook.

import pandas as pd
import numpy as np
import phasing as ph
import sh
import random
import os
import argparse

parser = argparse.ArgumentParser(description='Script to perform a shapeIt run')

parser.add_argument('input', help='input file root for .bed .bim .fam')
parser.add_argument('output', help='output directory')

#parser.add_argument('-Q', '--gq', dest='gq_threshold', action='store',
#                   help='Genotype quality filter to apply')

#parser.add_argument('-P', '--pedigree', dest='pedigree', action='store',
#                    help='Pedigree table for calculation of mendel errors')
args = parser.parse_args()

for ext in ('.bed', '.bim', '.fam'):
    assert os.path.isfile(args.input + ext)

# length of 3L
len3l = 41963435
segment_len = 5*(10**6)
overlap = 5*(10**5)

parameters = ['-B', args.input,
              '--thread', 4,
              '--no-ped',
              '--window', 0.5,
              '--rho', 0.0001,
              'effective-size', 500]

start_pos = 1
split_runs = dict()
while start_pos < len3l:
    stop_pos = start_pos + segment_len
    print(start_pos, stop_pos)

    shape_it = ph.algorithms.ShapeIt(parameters + ['--input-from', start_pos,
                                                   '--input-to', stop_pos],
                                     args.output)

    split_runs[shape_it.run_id] = shape_it
    shape_it.run('-l', 'h_vmem=8G', '-pe', 'simple_pe', '4')

    start_pos = stop_pos - overlap

# when all of these are finished we want to run the ligate haplotypes script

# and then maybe duoHMM

#du = ph.algorithms.DuoHMM(parameters=[], run_id=r.run_id,
                          #outdir=r.outdir,
                          #executable='/home/njh/exec/bin/duohmm')


#du.run('-l', 'h_vmem=4G', '-hold_jid', r.run_id)