#! /usr/bin/python
# simple script that takes an hdf5 file to create a PED file.
# note hdf5 files representing genotype data can be created using npyvcf

# input:
# - hdf5: including a calldata genotype field, representing the unphased genotypes.
# - tfam, a family file in transverse format as specified by: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr
# - region: a chr/region consideration. This also defines the name of the output dir.
# - output directory

import getopt
import sys
import pandas as pd
#import phasing.hdf5_2_tped

# Main exists to grab args. Calls create_tped_from_hdf5
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

    tfam = pd.read_csv(tfam_file, sep = " ", header = None)
    tfam = tfam.rename(columns = { 0:'Family_ID', 
                                   1:'Individual_ID', 
                                   2:'Paternal_ID', 
                                   3:'Maternal_ID', 
                                   4:'Sex', 
                                   5:'Phenotype'})

    #phasing.hdf5_2_tped.create_tped_from_hdf5(h5file, out_dir, region,
    # tfam.Individual_ID)
    
if __name__ == "__main__":
    main(sys.argv[1:])
