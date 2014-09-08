# Collection of functions for conversion of hdf5 variants to tped

# external
import numpy as np
import pandas as pd
import h5py
import os

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

def create_tped_from_hdf5(h5file, tfam_file, out_dir = None, region = None):

    # make checks
    assert os.path.isfile(h5file) 
    assert os.path.isfile(tfam_file) 

    contig, posrange = region.split(':')
    start, stop = [int(x) for x in posrange.split('-')]

    # load hdf5 file
    fh_hdf5 = h5py.File(h5file, mode = 'r')

    # load tfam_file_file
    h5_samples = fh_hdf5[contig]['samples'][:].tolist()

    tfam = pd.read_csv(tfam_file, sep = " ", header = None)
    tfam = tfam.rename(columns = { 0:'Family_ID', 
                                   1:'Individual_ID', 
                                   2:'Paternal_ID', 
                                   3:'Maternal_ID', 
                                   4:'Sex', 
                                   5:'Phenotype'})
    
    # read positions
    positions = fh_hdf5[contig]['variants']['POS'][:]
    variant_bool = np.array((positions >= start) & (positions <= stop))

    positions = np.compress(variant_bool, positions, axis=0)

    # read alt/ref
    ref = fh_hdf5[contig]['variants']['REF'][variant_bool]
    alt = fh_hdf5[contig]['variants']['ALT'][variant_bool]

    #include_samples = [s in tfam.Individual_ID for s in h5_samples]
    call_data = fh_hdf5[contig]['calldata']['genotype'][variant_bool,:,:]

    # determine which samples we care about: integer vector
    include_samples = [h5_samples.index(s) for s in tfam.Individual_ID]

    tped_file = os.path.join(out_dir, "_".join([ str(s) for s in [contig, start, stop]])+'.tped')
    fh_tped_file = open(tped_file, 'w');

    # loop through variants
    for idx in xrange(call_data.shape[0]):
        out_string = get_tped_row(call_data[idx, include_samples, :], 
                                  ref[idx], 
                                  alt[idx], 
                                  positions[idx], 
                                  contig)

        fh_tped_file.write(out_string + "\n");
    
    return tped_file
