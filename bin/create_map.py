__author__ = 'Nicholas Harding'

import numpy as np
import argparse
import time

start_time = time.time()
chunk_size = 2e5
parser = argparse.ArgumentParser(
    description='Tool to select the set of variants to phase from the'
                'truth set (from hdf5 file)')

parser.add_argument('output', help='output file')

# file: {OUT}/q{pacut}.q{prcut}.p{prop}.t{thin}/{FILESTEM}.h5
parser.add_argument('--geneticlength', '-M', action='store',
                    default=1.2, dest='geneticlength',
                    help='Length of chromosome in M')

parser.add_argument('--length', '-L', action='store',
                    default=41963435, dest='length',
                    help='Length of chromosome in bp')

parser.add_argument('--chr', '-C', action='store', default=None, dest='contig',
                    help='Which contig to evaluate.')

# offset is number by which we do not observe any variants here 10k
offset = 1e4

# to do: add option to only filter individual crosses.
args = parser.parse_args()
window_size = args.length / 500
positions = np.arange(offset, args.length + window_size, window_size).astype(
    "int")

with open(args.output, 'w') as f:
    f.write("\t".join(["pposition", "rrate", "gposition"]) + "\n")
    for pos in positions:
        gen_dist = (float(pos)/args.length) * args.geneticlength
        recomb_r = (args.geneticlength * 100) / (args.length/1e6)
        f.write("\t".join(str(s) for s in [pos, recomb_r, gen_dist]) + "\n")