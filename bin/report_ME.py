#! /usr/bin/python

__author__ = 'Nicholas Harding'

import h5py
import numpy as np
import anhima
import phasing as ph
import os
import argparse
import time


def me(pa, pr):
    return anhima.ped.diploid_mendelian_error(pa, pr).any(axis=1)


def li(pa, pr):
    return ph.utils.get_error_likelihood(pa, pr) < 0


start_time = time.time()
chunk_size = 2e5
parser = argparse.ArgumentParser(
    description='Tool to select the set of variants to phase from the'
                'truth set (from hdf5 file)')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('pedigree', help='path to load pedigree file')
parser.add_argument('outdir', help='output directory')

# arguments
# POSITONAL: input file, output directory
# ARGUMENTS:
# - pedigree

# file: {OUT}/q{pacut}.q{prcut}.p{prop}.t{thin}/{FILESTEM}.h5
parser.add_argument('--cross', '-P', action='append',
                    default=None, dest='cross',
                    help='Which pedigree to evaluate. Assumes all otherwise')

parser.add_argument('--chr', '-C', action='store', default=None, dest='contig',
                    help='Which contig to evaluate.')

parser.add_argument('--likelihood', action='store_const', default=me,
                    dest='bad_sites', const=li,
                    help='Which method by which to find bad sites.')

# to do: add option to only filter individual crosses.
args = parser.parse_args()
for arg, value in sorted(vars(args).items()):
    print arg, ":", str(value)

data = h5py.File(args.input)
pedigree, ped_table = ph.utils.read_pedigree_table(args.pedigree)

if args.cross is None:
    args.cross = pedigree.keys()

available_samples = data[args.contig]['samples'][:]

for x in args.cross:
    bad_sites = list()
    print 'Processing cross:', x
    assert np.in1d(pedigree[x]['parent'], available_samples).all(), \
        "not all parents are available for cross " + x
    assert np.in1d(pedigree[x]['progeny'], available_samples).all(), \
        "not all progeny are available for cross " + x

    fn = os.path.join(args.outdir, "_".join([x, args.contig,
                                             args.bad_sites.__name__,
                                             'badsites.npz']))
    k = args.contig
    max_position = data[k]['variants']['POS'][:].max()
    chunks = np.arange(1, max_position + chunk_size, chunk_size)
    assert chunks.max() > max_position

    #x_samples = pedigree[x]['parent'] + pedigree[x]['progeny']

    for start, stop in zip(chunks[:-1], chunks[1:]):
        var, par_calldata = anhima.h5.load_region(data, k, start, stop, ['POS'],
                                                  ['genotype'],
                                                  samples=pedigree[x]['parent'])

        _, pro_calldata = anhima.h5.load_region(data, k, start, stop, None,
                                                ['genotype'],
                                                samples=pedigree[x]['progeny'])

        errors = args.bad_sites(par_calldata['genotype'],
                                pro_calldata['genotype'])

        bad_sites.append(np.compress(errors, var['POS']))

    bad_sites = np.concatenate(bad_sites)
    print '{0}/{1} sites with errors found'.format(np.size(bad_sites),
                                                   max_position)
    np.savez_compressed(os.path.join(args.outdir, fn), positions=bad_sites)

print("--- Completed: {0} seconds ---".format(time.time() - start_time))