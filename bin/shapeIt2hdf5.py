__author__ = 'Nicholas Harding'

import argparse
from tables import *
import numpy as np
import gzip
from itertools import islice
import os

parser = argparse.ArgumentParser(description='Tool to convert shapeIt output '
                                             'to anhima readable')

parser.add_argument('haplotypes', help='haplotypes file')
parser.add_argument('samples', help='samples file')

parser.add_argument('-c', '--chunksize', dest='chunk_size', action='store',
                    type=int, default=100000,
                    help='Number of lines to read at once from file')

parser.add_argument('--chr', dest='chr', action='store',
                    help='contig label to use as key at root of hdf5')
parser.add_argument('-o', '--out', dest='outfile', action='store',
                    help='path to write hdf5 file')
args = parser.parse_args()

assert not os.path.isfile(args.outfile)
h5file = openFile(args.outfile, mode="w")

# Get the HDF5 root group
root = h5file.root

# Create the groups
chrom = h5file.create_group(root, args.chr)

grp_calldata = h5file.create_group(chrom, "calldata")
grp_variants = h5file.create_group(chrom, "variants")

position = h5file.create_earray(grp_variants, name='POS',
                                atom=IntAtom(itemsize=4),
                                expectedrows=1000000, shape=(0, ),
                                filters=Filters(complevel=9, complib='zlib'))

identif = h5file.create_earray(grp_variants, name='ID',
                               atom=StringAtom(itemsize=8),
                               expectedrows=1000000, shape=(0, ),
                               filters=Filters(complevel=9, complib='zlib'))

ref = h5file.create_earray(grp_variants, name='REF',
                           atom=StringAtom(itemsize=1),
                           expectedrows=1000000, shape=(0, ),
                           filters=Filters(complevel=9, complib='zlib'))
alt = h5file.create_earray(grp_variants, name='ALT',
                           atom=StringAtom(itemsize=1),
                           expectedrows=1000000, shape=(0, ),
                           filters=Filters(complevel=9, complib='zlib'))

f = gzip.open(args.haplotypes, 'rb')
# look at line to see how many samples...
n_haps = len(f.readline().split(' ')) - 5
f.seek(0, 0)  # reset reader

genotypes = h5file.create_earray(grp_calldata, name='genotype',
                                 atom=IntAtom(itemsize=1),
                                 expectedrows=1000000, shape=(0, n_haps/2, 2),
                                 filters=Filters(complevel=9, complib='zlib'))

while True:
    chunk = list(islice(f, args.chunk_size))
    if not chunk:
        break
    as_np = np.array([line.rstrip().split(' ') for line in chunk])
    position.append(as_np[:, 2].astype('int'))
    identif.append(as_np[:, 1])
    ref.append(as_np[:, 3])
    alt.append(as_np[:, 4])
    genotypes.append(as_np[:, 5:].astype('int').reshape((args.chunk_size,
                                                         n_haps/2, 2)))

h5file.close()