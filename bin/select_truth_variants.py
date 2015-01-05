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

# to do: add option to only filter individual crosses.
args = parser.parse_args()
for arg, value in sorted(vars(args).items()):
    print arg, ":", str(value)

data = h5py.File(args.input)

pedigree, ped_table = ph.utils.read_pedigree_table(args.pedigree)

samples = data['3L']['samples'][:].tolist()
position = data['3L']['variants']['POS'][:]
reference = data['3L']['variants']['REF'][:]
alternate = data['3L']['variants']['ALT'][:]

output_dir = os.path.join(args.outdir,
                          "q{0}_q{1}_p{2}_t{3}".format(args.parentGQ,
                                                       args.progenyGQ,
                                                       args.mincalled,
                                                       args.gap))
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

fn = os.path.join(output_dir, args.filestem + '.h5')
print "Will output to:", fn

for x in pedigree.keys():
    pedigree[x]['parent_idx'] = [samples.index(s)
                                 for s in pedigree[x]['parent']]
    pedigree[x]['progeny_idx'] = [samples.index(s)
                                  for s in pedigree[x]['progeny']]

# grab genotypes
genotypes = data['3L']['calldata']['genotype'][:]

# this could be sped up a little bit...
sites = list()
none_missing = list()

for x in pedigree.keys():

    parent_genotypes = genotypes[:, pedigree[x]['parent_idx']]
    parent_gqs = data['3L']['calldata']['GQ'][:, pedigree[x]['parent_idx']]

    heterozygotes = anhima.gt.is_het(parent_genotypes.any(axis=1))
    present_gts = anhima.gt.is_called(parent_genotypes).all(axis=1)

    meet_gq = np.all(parent_gqs >= args.parentGQ, axis=1)

    progeny_gq = data['3L']['calldata']['GQ'][:, pedigree[x]['progeny_idx']]
    prog_ok = np.mean(progeny_gq >= args.progenyGQ, axis=1) > 0.8
    sites.append(meet_gq & heterozygotes & prog_ok)
    none_missing.append(present_gts)

putative_sites = np.array(sites).T.any(axis=1)
none_missing_a = np.array(none_missing).T.all(axis=1)

keep = none_missing_a & putative_sites

print keep.sum(), '/', keep.size, 'sites meet all requirements of quality ' \
                                  'and are segregating'

# check for mendelian errors
mendelian_errors = list()
unlikely_genotypes = list()
for x in pedigree.keys():
    print 'Processing cross:', x
    genotypes_pa = genotypes[:, pedigree[x]['parent_idx']][keep]
    genotypes_pr = genotypes[:, pedigree[x]['progeny_idx']][keep]

    mendelian_errors.append(anhima.ped.diploid_mendelian_error(
        genotypes_pa, genotypes_pr).any(axis=1))

    unlikely_genotypes.append(
        ph.utils.get_error_likelihood(genotypes_pa, genotypes_pr) < 0)

mendelian_errors = np.vstack(mendelian_errors).T.any(axis=1)
unlikely_genotypes = np.vstack(unlikely_genotypes).T.any(axis=1)

keep_positions = np.compress(keep, position)

error_sites = np.compress(mendelian_errors | unlikely_genotypes, keep_positions)

# grab the positions of mendelian errors and unlikely
# genotypes to mask the testing data
np.savez_compressed(os.path.join(output_dir, args.filestem + '_badsites.npz'),
                    positions=error_sites)

print mendelian_errors.sum(), 'mendelian errors observed'
print unlikely_genotypes.sum(), 'unlikely parental genotypes observed'

# not wild about this code, but seems to work!
keep[keep] = np.invert(mendelian_errors | unlikely_genotypes)

# now do some thinning
keep_positions = np.compress(keep, position)

lp = keep_positions[0]
not_thinned = np.zeros(keep_positions.shape, dtype='bool')

for i, pos in enumerate(keep_positions[1:]):
    if pos >= lp + args.gap:
        not_thinned[i+1] = True
        lp = pos

# again, looks a little icky
keep[keep] = not_thinned
print "In total", keep.sum(), '/', keep.size, "sites retained"

output_h5 = h5py.File(fn, 'w-')

output_h5.create_group('/3L')
output_h5.create_group('/3L/variants')
output_h5.create_group('/3L/calldata')

output_h5.create_dataset(
    os.path.join('3L', 'variants', 'REF'),
    data=np.compress(keep, reference),
    chunks=(1000,),
    compression='gzip',
    compression_opts=1)

output_h5.create_dataset(
    os.path.join('3L', 'variants', 'ALT'),
    data=np.compress(keep, alternate),
    chunks=(1000,),
    compression='gzip',
    compression_opts=1)

output_h5.create_dataset(
    os.path.join('3L', 'variants', 'POS'),
    data=np.compress(keep, position),
    chunks=(1000,),
    compression='gzip',
    compression_opts=1)

output_h5.create_dataset(
    os.path.join('3L', 'samples'),
    data=np.array(samples))

output_h5.create_dataset(
    os.path.join('3L', 'calldata', 'genotype'),
    data=np.compress(keep, genotypes, axis=0),
    chunks=(1000, 10, 2),
    compression='gzip',
    compression_opts=1)
genotypes = None

genotype_qual = data['3L']['calldata']['GQ'][:]
output_h5.create_dataset(
    os.path.join('3L', 'calldata', 'GQ'),
    data=np.compress(keep, genotype_qual, axis=0),
    chunks=(1000, 10),
    compression='gzip',
    compression_opts=1)

output_h5.close()
print("--- Completed: {0} seconds ---".format(time.time() - start_time))