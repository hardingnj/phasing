__author__ = 'Nicholas Harding'

# This is a wrapper that combines the data from to datasets
# rewritten to use scikit-allele.
import argparse
from tables import Filters, openFile, IntAtom, StringAtom
import h5py
import allel
import numpy as np
from os.path import isfile
import time

start_time = time.time()
chunk_size = 200000

parser = argparse.ArgumentParser(
    description='Tool to produce a combined h5 file from 2 separate hdf5 files')

parser.add_argument('input1', help='input hdf5 file')
parser.add_argument('input2', help='input hdf5 file')
parser.add_argument('output', help='output file stem')

parser.add_argument('--contig', '-C', action='store',
                    dest='contig', required=True, type=str)
parser.add_argument('--samplesA', '-sA', action='store',
                    dest='samplesA', nargs="+")
parser.add_argument('--samplesB', '-sB', action='store',
                    dest='samplesB', nargs="+")

args = parser.parse_args()

print("--- Starting ---")
print(args)


def combine_sample_genotypes(path_x, path_y, output_path,
                             contig, samples_x=None, samples_y=None):

    if isfile(output_path):
        raise FileExistsError("out path already exists")

    h5file = openFile(output_path, mode="w")

    # load 1st hdf5
    fh_a = h5py.File(path_x)
    fh_b = h5py.File(path_y)

    # load genotypes
    ga = allel.GenotypeCArray.from_hdf5(fh_a[contig]["calldata"]["genotype"])
    gb = allel.GenotypeCArray.from_hdf5(fh_b[contig]["calldata"]["genotype"])

    alleles = ga.count_alleles()
    biallelic = np.array(alleles.max_allele() < 2)

    # load positions
    pos_a = fh_a[contig]["variants"]["POS"][:]
    pos_b = fh_b[contig]["variants"]["POS"][:]

    # filter out non-biallelic sites:
    ga = ga.compress(biallelic, axis=0)
    pos = np.compress(biallelic, pos_a, axis=0)
    ref = np.compress(biallelic, fh_a[contig]["variants"]["REF"][:], axis=0)
    alt = np.compress(biallelic, fh_a[contig]["variants"]["ALT"][:], axis=0)

    assert np.array_equal(pos, pos_b)

    # samples
    samplesa = fh_a[contig]["samples"][:]
    samplesb = fh_b[contig]["samples"][:]

    if samples_y:
        l = [s.decode() for s in samplesb]
        idx = [l.index(s) for s in samples_y]
        gb = gb.take(idx, axis=1)
        samplesb = samples_y

    if samples_x:
        l = [s.decode() for s in samplesa]
        idx = [l.index(s) for s in samples_x]
        ga = ga.take(idx, axis=1)
        samplesa = samples_x

    root = h5file.root

    # Create the groups
    chrom = h5file.create_group(root, contig)
    grp_calldata = h5file.create_group(chrom, "calldata")
    grp_variants = h5file.create_group(chrom, "variants")

    # create objects

    filters = Filters(complevel=1, complib='zlib')
    sample_names = np.concatenate([samplesa, samplesb]).astype("|S10")
    h5file.create_array(chrom, 'samples', sample_names)

    number_sites = ga.shape[0]

    position = h5file.create_earray(grp_variants, name='POS',
                                    atom=IntAtom(itemsize=4),
                                    expectedrows=number_sites, shape=(0, ),
                                    filters=filters)

    reference = h5file.create_earray(grp_variants, name='REF',
                                     atom=StringAtom(itemsize=1),
                                     expectedrows=number_sites, shape=(0, ),
                                     filters=filters)

    alternate = h5file.create_earray(grp_variants, name='ALT',
                                     atom=StringAtom(itemsize=1),
                                     expectedrows=number_sites, shape=(0, ),
                                     filters=filters)

    genotypes = h5file.create_earray(grp_calldata, name='genotype',
                                     atom=IntAtom(itemsize=1),
                                     expectedrows=number_sites,
                                     shape=(0, sample_names.size, 2),
                                     filters=filters)

    chunks = np.arange(0, number_sites, chunk_size)
    chunks[-1] = number_sites

    for start, stop in zip(chunks[:-1], chunks[1:]):

        gt = np.hstack([ga[start:stop], gb[start:stop]])
        genotypes.append(gt)

        position.append(pos[start:stop])
        reference.append(ref[start:stop])
        alternate.append(alt[start:stop])

    h5file.close()

combine_sample_genotypes(path_x=args.input1,
                         path_y=args.input2,
                         output_path=args.output,
                         contig=args.contig,
                         samples_x=args.samplesA,
                         samples_y=args.samplesB)

print("--- Completed: {0:.4} seconds ---".format(time.time() - start_time))