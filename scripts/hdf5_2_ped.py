#! /usr/bin/python
# simple script that takes an hdf5 file to create a PED file.
# note hdf5 files representing genotype data can be created using npyvcf

# input:
# - hdf5: including a calldata genotype field, representing the unphased genotypes.
# - tfam, a family file in transverse format as specified by: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr
# - region: a chr/region consideration. This also defines the name of the output dir.
# - output directory

import numpy as np
import pandas as pd
import h5py
import os
import sh
import sys
import getopt

def convert_gts_to_strings(genotypes, ref, alt):
    lu = { -1:'0', 0:ref, 1:alt }
    return np.apply_along_axis(lambda x: lu[x[0]] + ' ' + lu[x[1]], 1, genotypes).tolist()

def main(argv):

    h5file  = ''
    tfam_file = ''
    out_dir = None
    region  = None

    try:
        opts, args = getopt.getopt(
            argv,
            "i:t:o:r:",
            ["h5file=","tfam=", "outdir=", "region="])
    except getopt.GetoptError:
        print ' '.join([
            'hdf5_2_ped.py',
            '-i <inputfile (hdf5)>',
            '-t <tfam file>',
            '-o <outputdir>',
            '-r <region>'
        ])
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print ' '.join([
                'hdf5_2_ped.py',
                '-i <inputfile (hdf5)>',
                '-t <tfam file>',
                '-o <output directory>',
                '-r <region>'
            ])
            sys.exit()
        elif opt in ("-i", "--h5file"):
            h5file = arg
        elif opt in ("-o", "--outdir"):
            out_dir = arg
        elif opt in ("-t", "--tfam"):
            tfam_file = arg
        elif opt in ("-r", "--region"):
            region = arg

    # make checks
    assert os.path.isfile(h5file) 
    assert os.path.isfile(tfam_file) 

    contig, posrange = region.split(':')
    start, stop = [int(x) for x in posrange.split('-')]

    # load hdf5 file
    fh_hdf5 = h5py.File(h5file, mode = 'r')

    # load tfam_file_file
    h5_samples = fh_hdf5[contig]['samples'][:]

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

    # determine which samples we care about
    include_samples = [s in tfam.ID for s in h5_samples]

    call_data = fh_hdf5[contig]['calldata']['genotype'][variant_bool,:,:]

    tped_file = os.path.join(out_dir, "_".join([ str(s) for s in [contig, start, stop]])+'.tped')
    fh_tped_file = open(tped_file, 'w');

    # loop through variants
    for idx in xrange(call_data.shape[0]):
        str_gts = convert_gts_to_strings(call_data[idx, include_samples], ref[idx], alt[idx]);
        fh_tped_file.write("\t".join([contig, 'snp'+str(positions[idx]), '0', str(positions[idx])] + str_gts)+"\n")
    # end main
    
if __name__ == "__main__":
    main(sys.argv[1:])
