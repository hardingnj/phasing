__author__ = 'Nicholas Harding'

import phasing as ph
import os
import sh
import yaml
import argparse

parser = argparse.ArgumentParser(description='Script to run vcf creation '
                                             'pipeline')
parser.add_argument('config', help='setting file in yaml format')
args = parser.parse_args()
config = yaml.load(open(args.config, 'r'))


def create_command(raw=None, gq=None, out=None, ped=''):

    # summary
    # - filter hdf5 by GQ and optionally produce ME.
    # - convert to vcf format, skipping sites that failed to gt in all samples
    # - bgzip, tabix and remove h5 files.

    filter_command = config['pyenv'] + ' ' + config['bin_dir'] + \
        '/filter_hdf5.py {RAW} {OUT} -O -Q {GQ} {PED}'
    to_vcf_command = config['pyenv'] + ' ' + config['bin_dir'] + \
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
           rm_command.format(FILE=out+'.h5')]

    return cmd

# this python script supercedes the ipynb file vcf_to_plink
# This has several steps and could actually be done with a makefile. but as
# this requires splitting some files, I will do with file name checking.

# job dependencies not really handles currently. Just 2 independent script
# files.

# define raw hdf5 files.
hdf5_raw = {'eval': config['eval_h5'], 'truth': config['truth_h5']}

hdf5_stems = {k: os.path.splitext(os.path.basename(v))[0]
              for k, v in hdf5_raw.iteritems()}

# make the /vcf dir, script dir and the log dir in both eval/truth
assert os.path.isdir(config['outdir'])

dirs = dict()
for tpe in hdf5_raw.keys():
    dirs[tpe] = {d: os.path.join(config['outdir'], tpe, d)
                 for d in ('vcf', 'log', 'script')}
    [os.mkdir(d) for d in dirs[tpe].values() if not os.path.isdir(d)]

### TRUTH SET ##################################################################
ph.utils.create_sh_script(
    os.path.join(dirs['truth']['script'], 'prepare.sh'),
    create_command(raw=hdf5_raw['truth'], gq=config['gq_threshold'],
                   out=os.path.join(dirs['truth']['vcf'],
                                    hdf5_stems['truth']),
                   ped='-P ' + config['ped_file']))

sh.qsub('-l', config['truth_prep_mem'], '-N', 'prepare_truth', '-j', 'y',
        '-S', '/bin/bash', '-o', dirs['truth']['log'],
        os.path.join(dirs['truth']['script'], 'prepare.sh'))

### EVAL SET ###################################################################
ph.utils.create_sh_script(
    os.path.join(dirs['eval']['script'], 'prepare.sh'),
    create_command(raw=hdf5_raw['eval'], gq=config['gq_threshold'],
                   out=os.path.join(dirs['eval']['vcf'], hdf5_stems['eval'])))

sh.qsub('-l', config['eval_prep_mem'], '-N', 'prepare_eval', '-j', 'y',
        '-S', '/bin/bash', '-o', dirs['eval']['log'],
        os.path.join(dirs['eval']['script'], 'prepare.sh'))

# not add the ped file for the eval set.
# in both cases I want to calculate the PIRs. This can be left a bit.