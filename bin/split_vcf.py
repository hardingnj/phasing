import argparse
import phasing as ph
import os


def make_sample_file(sample_list, ped_dict, family_id, filepath):

    father, mother = ped_dict['parent']
    fh = open(filepath, 'w')
    fh.write('ID_1 ID_2 missing father mother sex plink_pheno' + "\n")
    fh.write('0 0 0 D D D B' + "\n")
    for i, s in enumerate(sample_list):
        if s in ped_dict['parent']:
            fh.write(" ".join([family_id + '_' + str(i+1), s, '0',
                               '0', '0', '0', '-9']) + "\n")
        elif s in ped_dict['progeny']:
            fh.write(" ".join([family_id + '_' + str(i+1), s, '0',
                               father, mother, '0', '-9']) + "\n")
        else:
            raise Exception('Sample not found in dict.')
    fh.close()

parser = argparse.ArgumentParser(description='Tool to filter a hdf5 file')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('output', help='output file stem')

parser.add_argument('-P', '--pedigree', dest='pedigree', action='store',
                    help='Pedigree table for calculation of mendel errors')
parser.add_argument('-B', '--binary', dest='binary', action='store',
                    default='vcfkeepsamples',
                    help='Path to vcfkeepsamples executable')
args = parser.parse_args()

pedigree, _ = ph.utils.read_pedigree_table(args.pedigree)

command_string = r"| perl -ne 'if ($_ =~ m/^#/ || $_ =~ m/0\/0|0\/1|1\/1/) " \
                 r"{ print $_;}' | bgzip -c >"

for k in pedigree.keys():
    samples = " ".join(pedigree[k]['parent'] + pedigree[k]['progeny'])
    print "Splitting " + k + ": " + samples
    cmd = " ".join([args.binary, args.input, samples, command_string,
                    args.output + '_' + k + '.vcf.gz'])
    os.system(cmd)

    make_sample_file(sample_list=samples.split(' '),
                     ped_dict=pedigree[k],
                     family_id=k,
                     filepath=args.output + '_' + k + '.sample')