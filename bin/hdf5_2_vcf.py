#! /usr/bin/python

__author__ = 'Nicholas Harding'

# This script takes an hdf5 file, reads the object in and produces a minimal
# vcf file. Really minimal, just the essential 8 columns and genotypes.

import argparse
import h5py
from itertools import izip
import anhima
import numpy as np
import phasing

chunk_size = 200000

parser = argparse.ArgumentParser(
    description='Tool to produce a vcf file from an hdf5 file')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('output', help='output file stem')

parser.add_argument('--keepmissing', '-M', action='store_true', default=False)
parser.add_argument('--cutoff', '-C', action='store', default=0.04,
                    dest='missingcutoff', type=float,
                    help='Maximum missing GTs tolerated in a sample')
parser.add_argument('--pedigree', '-P', action='store', dest='pedigree',
                    help='path to load pedigree file')

# to do: add option to only filter individual crosses.
args = parser.parse_args()

with h5py.File(args.input, mode='r') as h5_handle:
    with open(args.output + '.vcf', 'w') as f:

        lookup = {-1: './.', 0: '0/0', 1: '0/1', 2: '1/1'}

        f.write(r'##fileformat=VCFv4.1' + "\n")
        f.write(r'##FORMAT=<ID=GT,Number=1,Type=String,'
                r'Description="Genotype">' + "\n")
        f.write(r'##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count'
                r'in genotypes, for each ALT allele, in the same order as'
                r'listed">' + "\n")

        f.write(r'##contig=<ID=2L,length=49364325>' + "\n")
        f.write(r'##contig=<ID=2R,length=61545105>' + "\n")
        f.write(r'##contig=<ID=3L,length=41963435>' + "\n")
        f.write(r'##contig=<ID=3R,length=53200684>' + "\n")
        f.write(r'##contig=<ID=UNKN,length=42389979>' + "\n")
        f.write(r'##contig=<ID=X,length=24393108>' + "\n")
        f.write(r'##contig=<ID=Y_unplaced,length=237045>' + "\n")
        f.write(r'##reference=file:///data/anopheles/ag1000g/data/genome/AgamP3'
                r'/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.fa' + "\n")

        reqd = ('#CHROM', 'POS', 'ID', 'REF', 'ALT',
                'QUAL', 'FILTER', 'INFO', 'FORMAT')

        # rememeber to act on all 1st level keys!
        # does not support multiple chromosomes currently!
        # Actually should probably add to filter script...
        assert len(h5_handle.keys()) <= 1
        for k in h5_handle.keys():

            samples = tuple(h5_handle[k]['samples'][:].tolist())
            missing_rates = np.zeros(len(samples))

            for i, s in enumerate(samples):
                missing_genotypes = anhima.gt.is_missing(
                    h5_handle[k]['calldata']['genotype'][:, i].reshape(
                        (-1, 1, 2))).squeeze()
                consecutive_miss = phasing.utils.get_consecutive_true(
                    missing_genotypes)
                miss_rate_i = consecutive_miss/float(missing_genotypes.size)
                print "Missing rate of", s, ':', "{:.8f}".format(miss_rate_i), \
                    "({0}/{1})".format(i+1, len(samples))
                missing_rates[i] = miss_rate_i

            print "Rate max:", missing_rates.max()
            ok_samples = missing_rates < args.missingcutoff

            if np.any(~ok_samples):
                msg = "The following {0} samples are excluded as they have a " \
                      "a consecutive missing gt run of >= {1} of all calls:" \
                    .format(str(np.sum(~ok_samples)), str(args.missingcutoff))
                print msg

                for sa, rt in zip(
                        np.compress(~ok_samples, samples).tolist(),
                        np.compress(~ok_samples, missing_rates).tolist()):
                    print sa + ": " + str(rt)

                samples = tuple(np.compress(ok_samples, samples).tolist())
            else:
                print "All samples meet the missingness run threshold ({0})" \
                    .format(str(args.missingcutoff))

            if args.pedigree is not None:
                phasing.utils.create_samples_file(args.pedigree,
                                                  args.output + '.sample',
                                                  samples)

            f.write("\t".join(reqd + samples) + "\n")

            number_variants = h5_handle[k]['variants']['POS'][:].size
            chunks = np.arange(1, number_variants + chunk_size, chunk_size)
            assert chunks.max() > number_variants

            for start, stop in zip(chunks[:-1], chunks[1:]):
                sl = slice(start, stop)
                positions = h5_handle[k]['variants']['POS'][sl]
                reference = h5_handle[k]['variants']['REF'][sl]
                alternate = h5_handle[k]['variants']['ALT'][sl]
                genotypes = anhima.gt.as_012(
                    h5_handle[k]['calldata']['genotype'][sl])
                genotypes = np.compress(ok_samples, genotypes, axis=1)
                multiple_alts = alternate.ndim > 1

                for pos, ref, alt, gt in izip(positions, reference,
                                              alternate, genotypes):
                    filterstring = 'PASS'
                    # This line filters variants where ALL genotypes are missing
                    if not args.keepmissing and np.all(gt == -1):
                        continue

                    # alt may be an np array, with several entries.
                    if multiple_alts:
                        alt = ",".join(x for x in alt if x != '')

                    try:
                        line = "\t".join([k, str(pos), '.', ref, alt, '0', '.',
                                          '.', 'GT'] + [lookup[s] for s in gt])
                        f.write(line + "\n")

                    except TypeError:
                        print pos
                        print ref
                        print alt
                        print gt
                        raise TypeError("Some data wasn't of the correct type.")
