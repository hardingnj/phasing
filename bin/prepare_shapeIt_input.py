__author__ = 'Nicholas Harding'

# this python script supercedes the ipynb file vcf_to_plink that creates the
# ped files

# This has several steps and could actually be done with a makefile. but as
# this requires splitting some files, I will do with file name checking.

# job dependencies will be handles fine with hold_jid
# will try to split into mini scripts or even functions

# summary of pipeline:
# - input raw hdf5 file from VCF.

# 1. create TFAM file.

# 2. filter hdf5 of low GQ and mendelian errors, ie set as missing

# 3. convert hdf5 to tped files

# 4. convert tped files to binary format

# 5. merge into single AR bed/fim/fam file. Throw out variants that are missing
# in 10% or more of samples. --geno

# 6. Take the merged AR file, and convert it into ped format using plink.
# Then edit the plink file with perl to make all one family.
# also plit by family for easier merlin handling.

# Alternative pipeline:
# Use VCF that has been GQ filtered
# mask low qual genotypes using vcf tools
# create oxford/sample format using bcftools
# convert to ped/map using gtool
# convert to bed format
# create custom fam file with family information.
# filter out mendelian errors with plink.







base_dir = os.path.join('/data', 'anopheles', 'ag1000g', 'data', '1000g_09_13', 'evaluation',
                        'phasing', 'phase1.AR2', 'callsets', 'autosome', 'truth')
assert os.path.isdir(base_dir)

output_dir = os.path.join(base_dir, 'manipulated_formats')
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

print 'Output dir:\n', output_dir
pyexec = '/home/njh/pyenv/science/bin/python'

log_dir = os.path.join(output_dir, 'logs')
tped_dir = os.path.join(output_dir, 'tped')
script_dir = os.path.join(output_dir, 'scripts')
h5_dir = os.path.join(output_dir, 'h5')

if not os.path.isdir(log_dir):
    os.mkdir(log_dir)
if not os.path.isdir(tped_dir):
    os.mkdir(tped_dir)
if not os.path.isdir(script_dir):
    os.mkdir(script_dir)
if not os.path.isdir(h5_dir):
    os.mkdir(h5_dir)

import os
import sh
import re
import h5py
import pandas as pd
import phasing

# define hdf5 file position
h5input_dir = os.path.join('/data', 'anopheles', 'ag1000g', 'data', '1000g_09_13', 'evaluation',
                       'phasing', 'phase1.AR2', 'callsets', 'autosome', 'truth', 'vcf', 'env_f', 'h5')

hdf5_3L = os.path.join(h5input_dir, '3L_ag-cross.h5')
hdf5_3L_filtered = os.path.join(h5_dir, '3L_ag-cross.filtered.h5')

length_3L = 41963435
assert os.path.isfile(hdf5_3L)

# create tfam file: contains ID, pedigree, mother father, (optional sex and phenotype)
tfam_dir = os.path.join(output_dir, 'tfam')
if not os.path.isdir(tfam_dir):
    os.mkdir(tfam_dir)
tfam_file = os.path.join(tfam_dir, 'samples.tfam')

# Read in cross data
f_sample_ped = os.path.join(os.path.expanduser('~'), 'git', 'ag-crosses', 'meta', 'tbl_sample_ped.txt')
f_all_samples = os.path.join(os.path.expanduser('~'), 'git', 'ag-crosses', 'meta', 'all_samples.txt')

assert os.path.isfile(f_sample_ped) and os.path.isfile(f_all_samples)

ped_tbl = pd.read_csv(f_sample_ped, sep = "\t", index_col=1);
all_tbl = pd.read_csv(f_all_samples, sep = "\t", header = None)

all_tbl.columns = ['Individual_ID']
all_tbl.Individual_ID = [re.sub('_', '-', s) for s in all_tbl.Individual_ID]
all_tbl['pedigree'] = pd.Series();
all_tbl['Paternal_ID'] = -1;
all_tbl['Maternal_ID'] = -1;
all_tbl['Sex'] = 0
all_tbl['Phenotype'] = 0

for i, row in all_tbl.iterrows():
    try:
        match = ped_tbl.loc[row['Individual_ID']]
        cross = match['cross']
        all_tbl['pedigree'][i] = cross

        parents = ped_tbl[(ped_tbl.cross == cross) & (ped_tbl.role == 'parent')].index
        if match.name in parents:
            all_tbl['Paternal_ID'][i] = '0'
            all_tbl['Maternal_ID'][i] = '0'
        else:
            all_tbl['Paternal_ID'][i] = parents[0]
            all_tbl['Maternal_ID'][i] = parents[1]
    except:
        all_tbl['pedigree'][i] = 'none'
        all_tbl['Paternal_ID'][i] = '0'
        all_tbl['Maternal_ID'][i] = '0'

# Assign the family id. Basically {family}_{counter}
di = {}
l = []
for _, p in all_tbl.pedigree.iteritems():
    di[p] = di.get(p,0) + 1
    l.append(p + '_' + str(di[p]))
all_tbl['Family_ID'] = l

all_tbl[['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype']].to_csv(
    tfam_file,
    index = False,
    header = False,
    sep = " "
)
print 'tfam file written as:', tfam_file

# submit script
filter_script = os.path.join('/home', 'njh', 'git', 'phasing', 'bin', 'filter_hdf5.py')
filter_fn = os.path.join(script_dir, 'filter_hdf5_me.sh')
filter_log = os.path.join(log_dir, 'filter_hdf5_me.log')
filter_command = '{PYTHON} {SCRIPT} {IN} {OUT} -P {PED} -Q {GQ}'
filter_command = filter_command.format(PYTHON=pyexec,
                                       SCRIPT=filter_script,
                                       IN=hdf5_3L,
                                       OUT=hdf5_3L_filtered,
                                       PED=f_sample_ped,
                                       GQ=20)

phasing.utils.create_sh_script(filter_fn, [filter_command], hdf5_3L_filtered)
if not os.path.isfile(hdf5_3L_filtered + '.ok'):
    print sh.qsub('-l', 'h_vmem=8G','-N', 'hdf52_filter','-j', 'y', '-V',
                  '-o', filter_log, '-S', '/bin/bash',
                  filter_fn)

# ensure not used
hdf5_3L = None