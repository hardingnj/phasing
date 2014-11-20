import argparse
import h5py
import phasing as ph
import anhima
import numpy as np
import os


def copy_item(name, obj):
    group, dataset = os.path.split(name)
    if isinstance(obj, h5py.Dataset) and dataset != 'genotype':
        print 'copying ' + name + '...'
        filtered_h5.require_group(group)
        unfiltered_h5.copy(name, filtered_h5[group])

parser = argparse.ArgumentParser(description='Tool to filter a hdf5 file')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('output', help='output file stem')

parser.add_argument('-Q', '--gq', dest='gq_threshold', action='store',
                    type=int, help='Genotype quality filter to apply')

parser.add_argument('-P', '--pedigree', dest='pedigree', action='store',
                    help='Pedigree table for calculation of mendel errors')

parser.add_argument('-O', '--overwrite', dest='writemode', action='store_const',
                    const='w', default='w-',
                    help='Overwrites hdf5 file if already exists')
# to do: add option to only filter individual crosses.

args = parser.parse_args()

unfiltered_h5 = h5py.File(args.input, mode='r')
filtered_h5 = h5py.File(args.output + '.h5', mode=args.writemode)

# load ped file if defined
if args.pedigree is not None:
    pedigree, _ = ph.utils.read_pedigree_table(args.pedigree)

if args.pedigree is not None:
    fh = open(args.output + '_me.txt', args.writemode)

# remember to act on all 1st level keys!
for k in unfiltered_h5.keys():
    # copy groups to new object
    unfiltered_h5.visititems(copy_item)

    genotypes = unfiltered_h5[k]['calldata']['genotype'][:]
    positions = unfiltered_h5[k]['variants']['POS'][:]

    if args.gq_threshold is not None:
        genotype_qualities = unfiltered_h5[k]['calldata']['GQ'][:]
        mask = np.array(args.gq_threshold > genotype_qualities)
        assert mask.shape == genotype_qualities.shape
        genotypes[mask] = (-1, -1)

    samples = unfiltered_h5[k]['samples'][:].tolist()
    # now mask mendelian errors to an output file
    if args.pedigree is not None:
        # loop through crosses... and mark MEs
        mendelian_error = list()
        for c in pedigree:
            par_idx = np.array([samples.index(x) for x in pedigree[c][
                'parent']])
            prg_idx = np.array([samples.index(x) for x in pedigree[c][
                'progeny']])

            # mask bad parental genotypes
            parent_e = ph.utils.get_error_likelihood(genotypes[:, par_idx],
                                                     genotypes[:, prg_idx])
            genotypes[(parent_e < 0), par_idx] = (-1, -1)

            # mendelian errors:
            me = anhima.ped.diploid_mendelian_error(genotypes[:, par_idx],
                                                    genotypes[:, prg_idx])
            mendelian_error.append(np.any(me > 0, axis=1))




        variant_mask = np.vstack(mendelian_error).any(axis=0)

        mendel_err_positions = np.compress(variant_mask, positions)
        for x in mendel_err_positions:
            fh.write(k + "\t" + str(x) + "\n")
        fh.close()

    dset = filtered_h5.create_dataset(
        os.path.join(k, 'calldata', 'genotype'),
        data=genotypes,
        chunks=(1000, 10, 2),
        compression='gzip',
        compression_opts=1)