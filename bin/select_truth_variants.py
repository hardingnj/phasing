#! /usr/bin/python

__author__ = 'Nicholas Harding'


import h5py
import numpy as np
import anhima
import phasing as ph
import os
import argparse
import time

start_time = time.time()

parser = argparse.ArgumentParser(
    description='Tool to select the set of variants to phase from the'
                'truth set (from hdf5 file)')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('pedigree', help='path to load pedigree file')
parser.add_argument('outdir', help='output directory')

# arguments
# POSITONAL: input file, output directory
# ARGUMENTS:
# - filestem
# - %progeny not making filter
# - gq progeny cut off
# - gq parent cut off
# - thinning parameter

# file: {OUT}/q{pacut}.q{prcut}.p{prop}.t{thin}/{FILESTEM}.h5

parser.add_argument('--filestem', '-F', action='store',
                    default='segregating_sites')
parser.add_argument('--parentgq', '-Qp', action='store', default=40,
                    dest='parentGQ', type=int,
                    help='GQ filter applied to parents')
parser.add_argument('--progenygq', '-Qy', action='store', default=30,
                    dest='progenyGQ', type=int,
                    help='GQ filter applied to progeny')
parser.add_argument('--mincalled', '-M', action='store',
                    default=0.8, type=float, dest='mincalled',
                    help='minimum proportion of progeny genotypes that must be '
                         'non-missing')
parser.add_argument('--thinning', '-T', action='store', dest='gap',
                    default=100, type=int,
                    help='enforced gap between genotypes')

parser.add_argument('--baddir', '-B', action='store', dest='baddir',
                    default=os.getcwd(), type=str,
                    help='location of bad sites numpy files')

parser.add_argument('--cross', '-P', action='append',
                    default=None, dest='cross',
                    help='Which pedigree to evaluate. Assumes all otherwise')

parser.add_argument('--chr', '-C', action='store', default=None, dest='contig',
                    help='Which contig to evaluate.')


# to do: add option to only filter individual crosses.
args = parser.parse_args()
for arg, value in sorted(vars(args).items()):
    print arg, ":", str(value)

data = h5py.File(args.input)

pedigree, ped_table = ph.utils.read_pedigree_table(args.pedigree)
if args.cross is None:
    args.cross = pedigree.keys()

samples = data[args.contig]['samples'][:].tolist()
position = data[args.contig]['variants']['POS'][:]
reference = data[args.contig]['variants']['REF'][:]
alternate = data[args.contig]['variants']['ALT'][:]

output_dir = os.path.join(args.outdir,
                          "q{0}_q{1}_p{2}_t{3}".format(args.parentGQ,
                                                       args.progenyGQ,
                                                       args.mincalled,
                                                       args.gap))
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

fn = os.path.join(output_dir, args.filestem + '.h5')
output_h5 = h5py.File(fn, 'w-')
print "Will output to:", fn

for x in args.cross:
    pedigree[x]['parent_idx'] = [samples.index(s)
                                 for s in pedigree[x]['parent']]
    pedigree[x]['progeny_idx'] = [samples.index(s)
                                  for s in pedigree[x]['progeny']]

# grab genotypes
genotypes = data[args.contig]['calldata']['genotype'][:]

# this could be sped up a little bit...
hets = list()
min_quality = list()
bad_progeny = list()

for x in args.cross:

    parent_genotypes = genotypes[:, pedigree[x]['parent_idx']]
    parent_gqs = data[args.contig]['calldata']['GQ'][:, pedigree[x]['parent_idx']]

    # at least 1 parent must be a het
    heterozygotes = anhima.gt.is_het(parent_genotypes.any(axis=1))

    # both must be called
    present_gts = anhima.gt.is_called(parent_genotypes).all(axis=1)

    # both must be >= gq
    meet_gq = np.all(parent_gqs >= args.parentGQ, axis=1)

    progeny_gq = data[args.contig]['calldata']['GQ'][:, pedigree[x]['progeny_idx']]
    prog_bad = np.mean(progeny_gq >= args.progenyGQ, axis=1) < 0.8

    # append to lists
    hets.append(heterozygotes)                 # can be true in any 1
    min_quality.append(meet_gq & present_gts)  # must be true in all
    bad_progeny.append(prog_bad)               # must be false if parent is het

hets = np.array(hets).T
min_quality = np.array(min_quality).T
bad_progeny = np.array(bad_progeny).T

keep = hets.any(axis=1) & min_quality.all(axis=1) & ~np.any(hets &
                                                            bad_progeny, axis=1)

print '{0}/{1} sites meet all requirements of quality and are ' \
      'segregating'.format(np.sum(keep), np.size(keep))

# simply add all bad sites
bad_positions = list()
for x in args.cross:
    print 'Processing cross:', x
    for v in ['li', 'me']:
        fn = os.path.join(args.baddir,
                          "_".join([x, args.contig, v, 'badsites.npz']))
        bad_sites = np.load(fn)
        bad_positions.append(bad_sites['positions'])

bad_positions = np.unique(np.concatenate(bad_positions))

# whether pos is in bad pos list
error, _ = anhima.loc.locate_positions(position, bad_positions)
print '{0} sites blacklisted due to genotyping errors'.format(error.sum())
print '{0} sites removed with {1} remaining'.format(np.sum(error & keep),
                                                    np.sum(~error & keep))

keep = keep & ~error

# Thin the positions
keep_positions = np.compress(keep, position)

lp = keep_positions[0]    # bit tricky, won't allow passes until gets here
not_thinned = np.copy(keep)

for i in range(not_thinned.size - 1):
    if not not_thinned[i+1]:
        continue
    elif position[i+1] >= lp + args.gap:
        lp = position[i+1]
    else:
        not_thinned[i+1] = False

print "Thinning has removed {0} sites".format(np.sum(keep) - not_thinned.sum())
keep = not_thinned
print "In total {0}/{1} sites retained".format(np.sum(keep), np.size(keep))

output_h5.create_group('/3L')
output_h5.create_group('/3L/variants')
output_h5.create_group('/3L/calldata')

output_h5.create_dataset(
    os.path.join(args.contig, 'variants', 'REF'),
    data=np.compress(keep, reference),
    chunks=(1000,),
    compression='gzip',
    compression_opts=1)

output_h5.create_dataset(
    os.path.join(args.contig, 'variants', 'ALT'),
    data=np.compress(keep, alternate),
    chunks=(1000,),
    compression='gzip',
    compression_opts=1)

output_h5.create_dataset(
    os.path.join(args.contig, 'variants', 'POS'),
    data=np.compress(keep, position),
    chunks=(1000,),
    compression='gzip',
    compression_opts=1)

output_h5.create_dataset(
    os.path.join(args.contig, 'samples'),
    data=np.array(samples))

output_h5.create_dataset(
    os.path.join(args.contig, 'calldata', 'genotype'),
    data=np.compress(keep, genotypes, axis=0),
    chunks=(1000, 10, 2),
    compression='gzip',
    compression_opts=1)
genotypes = None

genotype_qual = data[args.contig]['calldata']['GQ'][:]
output_h5.create_dataset(
    os.path.join(args.contig, 'calldata', 'GQ'),
    data=np.compress(keep, genotype_qual, axis=0),
    chunks=(1000, 10),
    compression='gzip',
    compression_opts=1)

output_h5.close()
print("--- Completed: {0} seconds ---".format(time.time() - start_time))