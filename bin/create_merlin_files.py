__author__ = 'Nicholas Harding'


import os
import pandas as pd
import re
import argparse
import gzip


# This is inefficient- but does not involve reading huge amounts of data
# to memory
def get_ind_genotype(path, sample_index):
    convert = {'.': '0/0', './.': '0/0', '0/0': '1/1',
               '0/1': '1/2', '1/1': '2/2'}

    genotypes = list()
    fh = gzip.open(path, 'rb')
    for line_ in fh:
        if re.match('^#', line_):
            continue
        values = line_.rstrip().split("\t")
        try:
            genotypes.append(convert[values[9 + sample_index]])
        except KeyError:
            print genotypes
            print values
            exit(1)
    return "\t".join(genotypes)


parser = argparse.ArgumentParser(description='Tool to filter a hdf5 file')
parser.add_argument('vcf', help='input sample pedfile')
parser.add_argument('sample', help='File containing all samples')
parser.add_argument('output', help='output directory')

parser.add_argument('--recomb', '-r', action='store', dest='recomb',
                    default='110:40000000',
                    help='recomination rate in cM per chrom')
args = parser.parse_args()

recomb, length = [float(x) for x in args.recomb.split(':')]
cm_per_base = recomb / length

stem = os.path.split(args.vcf)[-1].split('.')[0]
print stem

# will produce 3 files: a ped, a dat and a map.

# create the dat and map files for merlin...
map_file = os.path.join(args.output, stem + '.map')
ped_file = os.path.join(args.output, stem + '.ped')
dat_file = os.path.join(args.output, stem + '.dat')

tbl_names = ('Family_ID', 'Individual_ID', 'missing', 'Paternal_ID',
             'Maternal_ID', 'Sex', 'Phenotype')
fam_tbl = pd.read_csv(args.sample, sep=" ", header=None,
                      names=tbl_names, skiprows=2)
fam_tbl['Family_ID'] = [re.sub('_\d+', '', x) for x in fam_tbl['Family_ID']]

vcf_fh = gzip.open(args.vcf, 'rb')
vcf_samples = None
for line in vcf_fh:
    if re.match('^##', line):
        continue
    elif re.match('^#CHR', line):
        vcf_samples = line.rstrip().split("\t")[9:]
        break

# dat file
dat_fh = open(dat_file, 'w')
dat_fh.write("\t".join(['T', 'phenotype']) + "\n")
map_fh = open(map_file, 'w')
map_fh.write("\t".join(['CHROMOSOME', 'MARKER', 'POSITION']) + "\n")
last_gp = 'X'

for line in vcf_fh:
    chrom, pos = line.split("\t")[:2]

    marker = 'id'+pos
    dat_fh.write("\t".join(['M', marker]) + "\n")

    gp = str(round(cm_per_base*int(pos), 6))
    try:
        assert gp != last_gp
    except AssertionError:
        print pos, last_gp
        Exception(AssertionError)
    map_fh.write("\t".join([chrom, marker, gp]) + "\n")

vcf_fh.close()
map_fh.close()
dat_fh.close()

# ped file
ped_fh = open(ped_file, 'w')
for i, s in enumerate(vcf_samples):
    match = fam_tbl[fam_tbl.Individual_ID == s]
    ped_fh.write("\t".join([str(x) for x in match.iloc[0].tolist()]) + "\t" +
                 get_ind_genotype(args.vcf, i) + "\n")
ped_fh.close()