__author__ = 'Nicholas Harding'

import phasing as ph
import os
import sh
import yaml
import argparse

parser = argparse.ArgumentParser(description='Script to run vcf creation '
                                             'pipeline')
parser.add_argument('config', help='setting file in yaml format')
parser.add_argument('-n', '--nosubmit', dest='nosubmit', action='store_const',
                    const=True, default=False,
                    help='If flag is passed script exits without running jobs')
args = parser.parse_args()
config = yaml.load(open(args.config, 'r'))


def split_vcfs_region(vcf_loc, vcf_stem, contig, start_pos, stop_pos, outdir):
    split_command = "vcftools --gzvcf {FILE} --out {OUT} --chr {CHR} " \
                    "--from-bp {START} --to-bp {STOP} --recode"
    bgzip_command = 'bgzip {OUT}.vcf'
    tabix_command = 'tabix -fp vcf {OUT}.vcf.gz'
    touch_command = 'touch {FILE}'

    vcf_out = os.path.join(outdir, vcf_stem + '.' +
                           "_".join([contig, str(start_pos), str(stop_pos)]))

    cmd = [split_command.format(FILE=os.path.join(vcf_loc,
                                                  vcf_stem + '.vcf.gz'),
                                OUT=vcf_out,
                                CHR=contig, START=start, STOP=stop),
           "mv {0} {1}".format(vcf_out + '.recode.vcf', vcf_out + '.vcf'),
           bgzip_command.format(OUT=vcf_out),
           tabix_command.format(OUT=vcf_out),
           touch_command.format(FILE=vcf_out + '.vcf.gz.ok')]

    return cmd


def create_vcf_cmd(raw=None, gq=None, out=None, ped=''):
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
    tabix_command = 'tabix -fp vcf {OUT}.vcf.gz'
    rm_command = 'rm {FILE}'

    cmd = [filter_command.format(RAW=raw, OUT=out, GQ=gq, PED=ped),
           touch_command.format(FILE=out+'.h5.ok'),
           to_vcf_command.format(H5=out+'.h5', OUT=out),
           bgzip_command.format(OUT=out),
           tabix_command.format(OUT=out),
           touch_command.format(FILE=out+'.vcf.gz.ok'),
           rm_command.format(FILE=out+'.h5')]

    return cmd


def mendelian_error_filter_cmd(vcf_stem, excl_sites):
    touch_command = 'touch {FILE}'
    bgzip_command = 'bgzip {OUT}.vcf'
    tabix_command = 'tabix -fp vcf {OUT}.vcf.gz'
    filter_me_command = "vcftools --gzvcf {VCF} --exclude-positions {FILE} " \
        "--recode --out {OUT}"

    cmd = [filter_me_command.format(VCF=vcf_stem + '.vcf.gz', FILE=excl_sites,
                                    OUT=vcf_stem),
           "mv {0} {1}".format(vcf_stem + '.recode.vcf',
                               vcf_stem + '.ME_filt.vcf'),
           bgzip_command.format(OUT=vcf_stem + '.ME_filt'),
           tabix_command.format(OUT=vcf_stem + '.ME_filt'),
           touch_command.format(FILE=vcf_stem + '.ME_filt.vcf.gz.ok')]

    return cmd

# this python script supercedes the ipynb file vcf_to_plink
# job dependencies not really handled currently. Just 2 independent script
# files.

# make the /vcf dir, script dir and the log dir in both eval/truth
assert os.path.isdir(config['outdir'])
#assert os.path.isfile(config['truth_h5'])
#assert os.path.isfile(config['eval_h5'])

eval_stem = os.path.splitext(os.path.basename(config['eval_h5']))[0] \
    + '.gq' + str(config['gq_threshold'])
truth_stem = os.path.splitext(os.path.basename(config['truth_h5']))[0] \
    + '.gq' + str(config['gq_threshold'])

eval_dirs = {d: os.path.join(config['outdir'], 'eval', d)
             for d in ('vcf', 'log', 'script', 'PIRs', 'vcf_regions')}

truth_dirs = {d: os.path.join(config['outdir'], 'truth', d)
              for d in ('vcf', 'log', 'script', 'PIRs', 'vcf_regions')}

[os.mkdir(d) for d in eval_dirs.values() + truth_dirs.values()
 if not os.path.isdir(d)]

truth_root_vcf = os.path.join(truth_dirs['vcf'], truth_stem)
eval_root_vcf = os.path.join(eval_dirs['vcf'], eval_stem)

make_vcf_name = str(config['gq_threshold']) + '_make_vcf.sh'
me_filter_name = str(config['gq_threshold']) + '_filter_me.sh'
split_vcf_name = str(config['gq_threshold']) + '_split_by_region.sh'

### 1 : CREATE VCF files
### 2 : FILTER BY ME SITES dependency on 1.
### 3 : SPLIT VCFS BY REGION
### 4 : RUN PIR JOBS

## TO DO: tidy up by adding a variable that is the stem, then add the ME_filt,
## and region to the end.

### 1 CREATE VCF FILES
# Truth
ph.utils.create_sh_script(
    os.path.join(truth_dirs['script'], make_vcf_name),
    create_vcf_cmd(raw=config['truth_h5'], gq=config['gq_threshold'],
                   out=truth_root_vcf, ped='-P ' + config['ped_file']))

# Evaluation
ph.utils.create_sh_script(
    os.path.join(eval_dirs['script'], make_vcf_name),
    create_vcf_cmd(raw=config['eval_h5'], gq=config['gq_threshold'],
                   out=eval_root_vcf))

### 2 FILTER MENDELIAN ERROR SITES
me_sites = os.path.join(truth_root_vcf + '_me.txt')

# Truth
ph.utils.create_sh_script(
    os.path.join(truth_dirs['script'], me_filter_name),
    mendelian_error_filter_cmd(
        vcf_stem=truth_root_vcf,
        excl_sites=me_sites))

# Evaluation
ph.utils.create_sh_script(
    os.path.join(eval_dirs['script'], me_filter_name),
    mendelian_error_filter_cmd(
        vcf_stem=eval_root_vcf,
        excl_sites=me_sites))

# in both cases I want to calculate the PIRs. This can be left a bit.

### 3 SPLIT RESULTING VCFs
window_regions = ph.utils.calc_regions(size=41963435, nbins=config['nbins'],
                                       overlap=None)

for i, region in enumerate(window_regions):
    start, stop = region
    task_id = str(i + 1)
    ph.utils.create_sh_script(
        os.path.join(truth_dirs['script'], task_id + '.' + split_vcf_name),
        split_vcfs_region(truth_dirs['vcf'], truth_stem + '.ME_filt', '3L',
                          start, stop, truth_dirs['vcf_regions']))
    ph.utils.create_sh_script(
        os.path.join(eval_dirs['script'], task_id + '.' + split_vcf_name),
        split_vcfs_region(eval_dirs['vcf'], eval_stem + '.ME_filt', '3L',
                          start, stop, eval_dirs['vcf_regions']))

### 4 SUBMIT PIR JOBS
# END
if args.nosubmit:
    exit(0)

### JOB SUBMISSION
# PHASE 1
if not os.path.isfile(truth_root_vcf + '.vcf.gz.ok'):
    sh.qsub('-l', config['truth_prep_mem'], '-N', 'prepare_truth',
            '-j', 'y', '-S', '/bin/bash', '-o', truth_dirs['log'],
            os.path.join(truth_dirs['script'], make_vcf_name))

if not os.path.isfile(eval_root_vcf + '.vcf.gz.ok'):
    sh.qsub('-l', config['eval_prep_mem'], '-N', 'prepare_eval', '-j', 'y',
            '-S', '/bin/bash', '-o', eval_dirs['log'],
            os.path.join(eval_dirs['script'], make_vcf_name))

# PHASE 2
if not os.path.isfile(truth_root_vcf + '.ME_filt.vcf.gz.ok'):
    sh.qsub('-l', 'h_vmem=4G', '-hold_jid', 'prepare_truth',
            '-N', 'me_filter_truth', '-j', 'y', '-S', '/bin/bash',
            '-o', truth_dirs['log'],
            os.path.join(truth_dirs['script'], me_filter_name))

if not os.path.isfile(truth_root_vcf + '.ME_filt.vcf.gz.ok'):
    sh.qsub('-l', 'h_vmem=4G', '-hold_jid', 'prepare_eval',
            '-N', 'me_filter_eval', '-j', 'y', '-S', '/bin/bash',
            '-o', truth_dirs['log'],
            os.path.join(eval_dirs['script'], me_filter_name))

# PHASE 3
sh.qsub('-l', 'h_vmem=2G', '-hold_jid', 'me_filter_truth',
        '-N', 'truth_split', '-j', 'y', '-S', '/bin/bash', '-b', 'y',
        '-o', truth_dirs['log'], '-t', '1-'+str(config['nbins']),
        os.path.join(truth_dirs['script'],
                     r"\\$SGE_TASK_ID." + split_vcf_name))
