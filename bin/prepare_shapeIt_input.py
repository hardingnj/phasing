__author__ = 'Nicholas Harding'

import phasing as ph
import os
import sh


def create_command(raw=None, gq=None, out=None, ped=''):

    # summary
    # - filter hdf5 by GQ and optionally produce ME.
    # - convert to vcf format, skipping sites that failed to gt in all samples
    # - bgzip, tabix and remove h5 files.

    filter_command = pyenv + ' ' + bin_dir + \
        '/filter_hdf5.py {RAW} {OUT} -O -Q ${GQ} {PED}'
    to_vcf_command = pyenv + ' ' + bin_dir + \
        '/hdf5_2_vcf.py {H5} {OUT}'
    touch_command = 'touch {FILE}'
    bgzip_command = 'bgzip {OUT}.vcf'
    tabix_command = 'tabix -p vcf {OUT}.vcf.gz'
    validate_command = 'vcf-validator {OUT}.vcf.gz'
    rm_command = 'rm {FILE}'

    cmd = [filter_command.format(RAW=raw, OUT=out, GQ=gq, PED=ped),
           touch_command.format(FILE=out+'.h5.ok'),
           to_vcf_command.format(H5=out+'.h5', OUT=out),
           bgzip_command.format(OUT=out),
           tabix_command.format(OUT=out),
           validate_command.format(OUT=out),
           rm_command.format(FILE=out+'h5')]

    return cmd

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


# set some parameters. Probably will set these a config file.
gq_threshold = 30
pyenv = '/home/njh/pyenv/science/bin/python'
bin_dir = '/home/njh/git/phasing/bin'
ped_file = '/home/njh/git/ag-crosses/meta/tbl_sample_ped.txt'

# raw definitions
base_outdir = os.path.join('/data', 'anopheles', 'ag1000g', 'data',
                           '1000g_09_13', 'evaluation', 'phasing', 'phase1.AR2',
                           'callsets', 'autosome')
assert os.path.isdir(base_outdir)

# define raw hdf5 files.
hdf5_raw = dict()
hdf5_raw['eval'] = os.path.join('/data/anopheles/ag1000g/data/1000g_09_13/',
                                'release/AR2/ag1000g.phase1.AR2.3L.PASS.h5')
hdf5_raw['truth'] = os.path.join('/data/anopheles/ag-crosses/data/release/',
                                 '0.1.GATK.PHASING.1000AG.AR2/h5/',
                                 '3L_ag-cross.h5')

# handle the truth (can I do this??)
# make the /vcf dir, script dir and the log dir in both eval/truth
dirs = dict()
for tpe in hdf5_raw.keys():
    dirs[tpe] = {d: os.path.join(base_outdir, tpe, d)
                 for d in ('vcf', 'log', 'script')}
    [os.mkdir(d) for d in dirs[tpe].values() if not os.path.isdir(d)]

# in both cases I want to run the filter script, then the 2 vcf script, but
ph.utils.create_sh_script(os.path.join(dirs['truth']['script'], 'prepare.sh'),
                          create_command(raw=hdf5_raw['truth'],
                                         gq=gq_threshold,
                                         out=dirs['truth']['vcf'],
                                         ped='-P ' + ped_file))

sh.qsub('-l', 'h_vmem=16G', '-N', 'prepare_truth', '-j', 'y',
        '-S', '/bin/bash', '-o', dirs['truth']['log'],
        os.path.join(dirs['truth']['script'], 'prepare.sh'))

# not add the ped file for the eval set.

# in both cases I want to calculate the PIRs. This can be left a bit.

