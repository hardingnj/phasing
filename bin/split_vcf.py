import argparse
import phasing as ph
import os

parser = argparse.ArgumentParser(description='Tool to filter a hdf5 file')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('output', help='output file stem')

parser.add_argument('-P', '--pedigree', dest='pedigree', action='store',
                    help='Pedigree table for calculation of mendel errors')

args = parser.parse_args()

pedigree, _ = ph.utils.read_pedigree_table(args.pedigree)

command_string = r"vcfkeepsamples {FILE} {SAMPLES} | perl -ne 'if ($_ =~ m/^#/ || $_ =~ m/0\/0|0\/1|1\/1/) { print $_;}' | bgzip -c >  {OUTPUT}.vcf.gz"

for k in pedigree.keys():
    samples = " ".join([pedigree[k]['parent'] + pedigree[k]['progeny'])
	cmd = command_string.format(FILE=args.input,
                                SAMPLES=samples,
                                OUTPUT=args.output + '_' + k + '.vcf.gz')
    os.system(cmd)
