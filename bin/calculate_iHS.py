__author__ = 'Nicholas Harding'

import h5py
import numpy as np
import argparse
import time
import allel
import pandas as pd
from os.path import join
# walk along the snps, until we see that the EHH has dropped to 0.05.
# So We need the EHH at every point.

from itertools import combinations


def calculate_EHH(block):
    # takes 1 block. Calculates the pairwise mean that they are identical.
    return np.mean([np.array_equal(block.T[i], block.T[j])
                    for i, j in combinations(range(block.shape[1]), 2)])


def chrom_walk(hptyp, start_pos, end_pos, end_when=0.0, buffer=False):
    invariant = np.all(hptyp == 0, axis=1) | np.all(hptyp == 1, axis=1)

    di = (end_pos - start_pos)//abs(end_pos - start_pos)
    ehh_list = [1.0]
    for x in np.arange(start_pos + di, end_pos, di):
        # no point calculating if invariant, same as last entry.
        if invariant[x]:
            ehh_list.append(ehh_list[-1])
            continue
        else:
            ehh = calculate_EHH(hptyp[start_pos:x:di])
            ehh_list.append(ehh)
            if ehh <= end_when:
                break

    ehh_out = np.array(ehh_list)
    if buffer:
        # we fill the rest of the array with zeros
        short = abs(end_pos - start_pos) - ehh_out.size
        ehh_out = np.append(ehh_out, np.zeros(short))

    return ehh_out


def calculate_iHS(locations, pos, haplotypes, fn):
    # takes an array of ps, the position array, and the gt array.

    assert haplotypes.ndim == 2
    assert isinstance(haplotypes, np.ndarray), "haps must be a np array"
    assert pos.ndim == 1
    print(haplotypes.shape)

    # for each "p"
    with open(fn, "w") as wp:

        for p in locations:
            core_anc = np.where(haplotypes[p] == 0)[0]
            core_der = np.where(haplotypes[p] == 1)[0]
            der_freq = np.mean(haplotypes[p] == 1)

            # divides it into coreA and B. idx
            core_anc_hap = haplotypes.take(indices=core_anc, axis=1)
            core_der_hap = haplotypes.take(indices=core_der, axis=1)

            # for coreB (derived) moves left and right calculating EHH until
            # hits < 0,05
            right_ehh = chrom_walk(core_der_hap, p, haplotypes.shape[0], 0.05)
            left_ehh = chrom_walk(core_der_hap, p, 0, 0.05)

            right_anc_ehh = chrom_walk(core_anc_hap, p, p + right_ehh.size,
                                       buffer=True)
            left_anc_ehh = chrom_walk(core_anc_hap, p, p - left_ehh.size,
                                      buffer=True)

            # now calculate trapezoid
            # right

            r_positions = pos[p:p + right_ehh.size]
            r_der_int = np.trapz(y=right_ehh, x=r_positions)
            r_anc_int = np.trapz(y=right_anc_ehh, x=r_positions)

            l_positions = pos[p:p - left_ehh.size: -1]
            l_der_int = -np.trapz(y=left_ehh, x=l_positions)
            l_anc_int = -np.trapz(y=left_anc_ehh, x=l_positions)
            unadjusted_ihs = np.log((r_der_int + l_der_int) /
                                    (r_anc_int + l_anc_int))

            d = [str(s) for s in (pos[p], der_freq,
                 r_der_int + l_der_int,
                 r_anc_int + l_anc_int,
                 unadjusted_ihs)]

            print("\t".join(d), file=wp)

    return None

start_time = time.time()

parser = argparse.ArgumentParser(
    description='Tool to identify areas under selection according to iHS')

# files:
parser.add_argument('input', help='input hdf5 filestem of phased data')
parser.add_argument('output', help='output directory')

parser.add_argument('--meta', '-M', action='store', default=None,
                    dest='meta', required=True, help="filepath for sample info")

# data:
parser.add_argument('--population', '-P', action='store', default=None,
                    dest='population', help='Which population to evaluate.')

parser.add_argument('--chr', '-C', action='store', default=None,
                    dest='chrom', help='Which contig to analyse')

# parameters
parser.add_argument('--region', '-L', action='store', dest='region',
                    type=str, help='Start-Stop position')

args = parser.parse_args()

chrom = args.chrom
fh = h5py.File(args.input, "r")[chrom]
g = allel.GenotypeCArray.from_hdf5(fh['calldata']['genotype'])
fh_samples = [s.decode() for s in fh['samples'][:]]

dat = pd.read_csv(args.meta, sep="\t", index_col=1)

samples = dat.loc[dat.population == args.population].index.tolist()

sample_idx = [fh_samples.index(s) for s in samples]

gt = g.take(sample_idx, axis=1)

all_inds_hom = np.all(gt.is_hom_ref(), axis=1) | np.all(gt.is_hom_alt(), axis=1)

gt = gt.compress(~all_inds_hom)
positions = np.compress(~all_inds_hom, fh['variants']['POS'][:])
alleles = gt.count_alleles()

max_allele = alleles.max_allele()
biallelic = max_allele <= 1
#allelism = alleles.allelism()

af = np.array(alleles.to_frequencies())
crit = af.T[0:2].min(axis=0) > 0.05
#is_seg = allelism > 1

to_analyse = crit & biallelic
loci = np.where(to_analyse)[0]
pos_loci = np.take(positions, loci)

start, stop = [int(q) for q in args.region.split("-")]
print("REGION:", start, stop)

idx = np.searchsorted(pos_loci, [start, stop])
loci_to_analyse = loci[slice(*idx)]
print("Analysing:", np.take(positions, loci_to_analyse))
print("nSites:", loci_to_analyse.size)

haps = np.array(gt.to_haplotypes())

output_fn = join(args.output, "{chrom}_{pop}_{start}_{stop}_iHS.txt".format(
    chrom=chrom, pop=args.population, start=start, stop=stop))
calculate_iHS(loci_to_analyse, positions, haps, output_fn)

end_time = time.time()
duration = (end_time - start_time)//60

print("----------------------------------")
print("Took {0} minutes.".format(duration))
print("----------------------------------")

