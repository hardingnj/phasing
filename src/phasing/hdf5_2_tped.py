# Collection of functions for conversion of hdf5 variants to tped

# external
import numpy as np
import h5py
import os
import anhima.loc

# This function returns a list of genotypes with alleles separated by spaces. Missing values are encoded as zeroes.
# - genotypes: a numpy 2D array of N x 2 genotypes as encoded by anhima. ie 0/1/-1 is ref/alt/missing
# - ref      : reference base at this position
# - alt      : alternate base at this position

# - output   : string list of length 'N'. Alleles separated by a space

def convert_gts_to_strings(genotypes, ref, alt):

    lu = { -1:'0', 0:ref, 1:alt }
    return np.apply_along_axis(lambda x: lu[x[0]] + ' ' + lu[x[1]], 1, genotypes).tolist()

# Given a genotype array, which samples, a ref, alt, pos, and contig this function will return the correct tped line.

def get_tped_row(gt_data, reference, alternate, position, contig):

        str_gts = convert_gts_to_strings(gt_data, reference, alternate);
        return "\t".join([contig, 'snp'+str(position), '0', str(position)] + str_gts)


def create_tped(gtypes, reference, alternate, position, fn_tped, contig = '0'):

    fh_tped_file = open(fn_tped, 'w');

    # loop through variants
    for idx in xrange(gtypes.shape[0]):
        out_string = get_tped_row(gtypes[idx, :, :], 
                                  reference[idx], 
                                  alternate[idx], 
                                  position[idx], 
                                  contig)

        fh_tped_file.write(out_string + "\n");
    
    return None

def create_tped_from_hdf5(h5file, out_dir = None, region = None, samples = None):

    # make checks
    assert os.path.isfile(h5file) 

    contig, posrange = region.split(':')
    start, stop = [int(x) for x in posrange.split('-')]

    # load hdf5 file
    fh_hdf5 = h5py.File(h5file, mode = 'r')

    # determine samples that we will uses
    h5_samples = fh_hdf5[contig]['samples'][:].tolist()

    if samples is None:
        samples = h5_samples

    # determine which samples we care about: integer vector
    include_samples = [h5_samples.index(s) for s in samples]

    # read positions
    pos = fh_hdf5[contig]['variants']['POS'][:]

    position_slice = anhima.loc.locate_region(pos,
                                              start_position=start, 
                                              stop_position=stop)

    pos = pos[position_slice]

    # read alt/ref
    ref = fh_hdf5[contig]['variants']['REF'][position_slice]
    alt = fh_hdf5[contig]['variants']['ALT'][position_slice]

    call_data = fh_hdf5[contig]['calldata']['genotype'][position_slice,:,:]
    tped_file = os.path.join(out_dir, "_".join([ str(s) for s in [contig, start, stop]])+'.tped')

    create_tped(call_data[:,include_samples,:], ref, alt, pos, tped_file, contig)

    return tped_file
