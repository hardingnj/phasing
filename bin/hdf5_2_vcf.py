__author__ = 'Nicholas Harding'

# This script takes an hdf5 file, reads the object in and produces a minimal
# vcf file. Really minimal, just the essential 8 columns and genotypes.

import argparse
import h5py
from itertools import izip
import anhima
import numpy as np

parser = argparse.ArgumentParser(
    description='Tool to produce a vcf file from an hdf5 file')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('output', help='output file stem')

parser.add_argument('--keepmissing', '-m', action='store_true', default=False)

# to do: add option to only filter individual crosses.
args = parser.parse_args()

h5_handle = h5py.File(args.input, mode='r')

lookup = {-1: './.', 0: '0/0', 1: '0/1', 2: '1/1'}

f = open(args.output + '.vcf', 'w')

f.write(r'##fileformat=VCFv4.1' + "\n")
f.write(r'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
f.write(r'##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in '
        r'genotypes, for each ALT allele, in the same order as listed">' +
        "\n")

f.write(r'##contig=<ID=2L,length=49364325>' + "\n")
f.write(r'##contig=<ID=2R,length=61545105>' + "\n")
f.write(r'##contig=<ID=3L,length=41963435>' + "\n")
f.write(r'##contig=<ID=3R,length=53200684>' + "\n")
f.write(r'##contig=<ID=UNKN,length=42389979>' + "\n")
f.write(r'##contig=<ID=X,length=24393108>' + "\n")
f.write(r'##contig=<ID=Y_unplaced,length=237045>' + "\n")
f.write(r'##reference=file:///data/anopheles/ag1000g/data/genome/AgamP3'
        r'/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa' + "\n")

reqd = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')

# rememeber to act on all 1st level keys!
for k in h5_handle.keys():

    genotypes = anhima.gt.as_012(h5_handle[k]['calldata']['genotype'][:])
    samples = tuple(h5_handle[k]['samples'][:].tolist())
    f.write("\t".join(reqd + samples) + "\n")

    positions = h5_handle[k]['variants']['POS'][:]
    reference = h5_handle[k]['variants']['REF'][:]
    alternate = h5_handle[k]['variants']['ALT'][:]

    for pos, ref, alt, gt in izip(positions, reference, alternate, genotypes):
        filterstring = 'PASS'

        # This line should filter out variants where ALL genotypes are missing
        if not args.keepmissing and np.all(gt == -1):
            continue

        try:
            line = "\t".join([k, str(pos), '.', ref, alt, '0', '.', '.',
                              'GT'] + [lookup[int(s)] for s in gt])
            f.write(line + "\n")
        except TypeError:
            print pos
            print ref
            print alt
            print gt
            exit(1)


f.close()
