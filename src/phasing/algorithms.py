__author__ = 'Nicholas Harding'

import re
import os
import tool
import pandas as pd
import numpy as np
import sh
import uuid
import tempfile
import utils
import subprocess
import yaml
# collection of methods for submission of shapeit jobs with various defaults

# behaviour:

# constructor will work out version
# print out handy yaml file with settings etc
# the class represents a run and keeps track of all the files etc.

# main function returns a string with the shapeit command

# also have a method that submits a qsub. This is generic

# methods that create the dir structure etc. Probably generic too

# data/shapeit/version/

# what would happen if I make duoHMM inherit from shape it:
# does it hold for multiple algs?
# better to pass


# needs
# would like shapeIt run to be aware that duoHMM has been run on it, so can
# parse results
# duo hmm needs to be aware of the parameters of shape it
# add duo hmm as a flag to shape It.
def parse_command(parameters):

    cl = re.compile('^-')
    last_was_key = False
    key = None
    parameters = [str(x) for x in parameters]

    command_dict = {}
    for value in parameters:
        value = str(value)
        if cl.match(value):
            key = value
            command_dict[key] = True
            last_was_key = True
        else:
            if last_was_key:
                command_dict[key] = value
                last_was_key = False
            else:
                command_dict[key] = command_dict[key] + ';' + value

    return command_dict


class ShapeIt():

    @staticmethod
    def _get_version(executable):

        p = subprocess.Popen(args=[executable, '--version'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, _ = p.communicate()

        p = re.compile("Version : (.+)\n")  # parentheses for capture groups
        m = p.search(out)
        if m:
            return m.group(1)
        else:
            print(out)
            raise Exception('Version not parsed.')

    def __init__(self, outdir, executable='shapeit',
                 ligate_exec='/home/njh/exec/ligateHAPLOTYPES/bin'
                             '/ligateHAPLOTYPES',
                 duohmm_exec='/home/njh/exec/bin/duohmm',
                 version=None):

        self.executable = executable
        self.ligate_bin = ligate_exec
        self.duohmm_bin = duohmm_exec
        if version is None:
            version = self._get_version(self.executable)
        self.version = version
        self.name = 'ShapeIt'
        self.checksum_file = None

        self.run_id = self.name + '_' + str(uuid.uuid4().get_hex().upper()[0:8])

        self.outdir = os.path.join(outdir, self.name, self.version, self.run_id)
        self.basedir = outdir
        self.si_job_list = list()
        self.ligate_script = None
        self.duohmm_script = None
        self.si_script = None
        self.dirs = {d: os.path.join(self.outdir, d) for d in ('log', 'script')}
        [sh.mkdir(d, '-p') for d in self.dirs.values() if not os.path.isdir(d)]

        self.haplotypes_f = os.path.join(self.outdir, self.run_id +
                                         '_final.haps.gz')
        self.phased_f = os.path.join(self.outdir, self.run_id +
                                     '_final.samples.gz')

        self.duohmm_haps = None
        self.duohmm_sample = None

        self.param_f = os.path.join(self.outdir, self.run_id + '.yaml')
        self.settings = {'run_id': self.run_id, 'version': self.version,
                         'executable': self.executable, 'outdir': self.outdir,
                         'basedir': outdir, 'haps': self.haplotypes_f,
                         'samples': self.phased_f}

    def setup_region_jobs(self, parameters, regions=None, vcf_file=None,
                          pirs=None, duohmm=False, sample_file=None):

        # basically, set up a shapeIT job for each region. This may be several
        # lines long.
        region_dir = os.path.join(self.outdir, 'region_outputs')
        os.mkdir(region_dir)

        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file

        if pirs is not None:
            parameters.insert(0, '-assemble')

        hap_files = list()
        for i, region in enumerate(regions):
            start, stop = [str(x) for x in region]
            haps = os.path.join(self.outdir, region_dir, str(i) + '_' +
                                self.run_id + '.haps.gz')
            hap_files.append(haps)
            samples = os.path.join(self.outdir, region_dir, str(i) + '_' +
                                   self.run_id + '.sample.gz')

            pir_string = ''
            if pirs is not None:
                pir_f = pirs.format(start, stop)
                assert os.path.isfile(pir_f)
                pir_string = " ".join(['--input-pir', pir_f])
                print pir_string

            files = ['--input-from', start, '--input-to', stop, pir_string,
                     '--output-max', haps, samples]
            print files

            cmd_shape_it = " ".join([self.executable] + parameters + files)

            script_name = os.path.join(self.dirs['script'], str(i) +
                                       '_shapeIt.sh')
            exit_check = "\n".join(["{",
                                    "if [ -f %s.ok ]; then",
                                    "\techo \"Already completed ok!\"",
                                    "\texit 0",
                                    "fi",
                                    "}"])
            utils.create_sh_script(
                filename=script_name,
                commands=['cd ' + region_dir, exit_check % haps,
                          cmd_shape_it],
                outfile=haps)

            # name, script, mem, dependency
            self.si_job_list.append(script_name)

        # set up ligateHaplotypes
        tmp = tempfile.NamedTemporaryFile(delete=False)
        gunzip = "gunzip -c {0} > {1}"

        self.ligate_script = os.path.join(self.dirs['script'], 'ligatehaps.sh')
        cmd_ligate = [self.ligate_bin, '--vcf', tmp.name, '--chunks',
                      " ".join(hap_files), '--output', self.haplotypes_f,
                      self.phased_f]

        utils.create_sh_script(self.ligate_script,
                               [gunzip.format(vcf_file, tmp.name),
                                'cd ' + self.outdir,
                                " ".join(cmd_ligate),
                                'rm ' + tmp.name],
                               self.haplotypes_f)

        if duohmm:
            tmp_duo = tempfile.NamedTemporaryFile(delete=False)
            self.duohmm_script = os.path.join(self.dirs['script'],
                                              'duohmm.sh')

            duohmm_root = os.path.join(self.outdir, self.run_id + '_duohmm')
            self.duohmm_haps = duohmm_root + '.haps'
            self.duohmm_sample = duohmm_root + '.sample'

            cmd_duohmm = [self.duohmm_bin, '-H', tmp_duo.name,
                          '-O', duohmm_root,
                          '-G', duohmm_root + '.GE.txt',
                          '-R', duohmm_root + '.RC.txt']

            utils.create_sh_script(self.duohmm_script,
                                   [gunzip.format(self.haplotypes_f,
                                                  tmp_duo.name+'.haps'),
                                   "cp {0} {1}".format(sample_file,
                                                       tmp_duo.name+'.sample'),
                                   'cd ' + self.outdir,
                                   " ".join(cmd_duohmm),
                                   "gzip {0} {1}".format(self.duohmm_haps,
                                                         self.duohmm_sample)],
                                   self.duohmm_haps + '.gz')

        self.settings['params'] = parse_command(parameters)

    def run_region(self, si_args, li_args, dm_args=None):

        qsub_parameters = ['-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.dirs['log']]
        # get checksum
        self.settings['checksum'] = utils.md5_for_file(self.checksum_file)
        yaml.dump(self.settings,
                  stream=open(self.param_f, 'w'))

        for job in self.si_job_list:
            print sh.qsub('-N', self.run_id, qsub_parameters, si_args, job)

        print sh.qsub('-hold_jid', self.run_id, '-N', 'ligate' + self.run_id,
                      qsub_parameters + li_args, self.ligate_script)

        if self.duohmm_script is not None:
            print sh.qsub('-hold_jid', 'ligate' + self.run_id,
                          '-N', 'duohmm' + self.run_id,
                          qsub_parameters + dm_args, self.duohmm_script)

    def setup_single_job(self, parameters, vcf_file, duohmm, sample_file=None):
        # basically, set up a shapeIT job running the whole file as one.
        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file

        cmd_shape_it = " ".join([self.executable] + parameters +
                                ['--output-max', self.haplotypes_f,
                                 self.phased_f])

        self.si_script = os.path.join(self.dirs['script'], 'shapeIt.sh')
        utils.create_sh_script(filename=self.si_script,
                               commands=['cd ' + self.outdir,
                                         cmd_shape_it],
                               outfile=self.haplotypes_f)

        if duohmm:
            gunzip = "gunzip -c {0} > {1}"
            tmp_duo = tempfile.NamedTemporaryFile(delete=False)
            self.duohmm_script = os.path.join(self.dirs['script'],
                                              'duohmm.sh')

            duohmm_root = os.path.join(self.outdir, self.run_id + '_duohmm')
            self.duohmm_haps = duohmm_root + '.haps'
            self.duohmm_sample = duohmm_root + '.sample'

            cmd_duohmm = [self.duohmm_bin, '-H', tmp_duo.name,
                          '-O', duohmm_root,
                          '-G', duohmm_root + '.GE.txt',
                          '-R', duohmm_root + '.RC.txt']

            utils.create_sh_script(self.duohmm_script,
                                   [gunzip.format(self.haplotypes_f,
                                                  tmp_duo.name+'.haps'),
                                   "cp {0} {1}".format(sample_file,
                                                       tmp_duo.name+'.sample'),
                                   'cd ' + self.outdir,
                                   " ".join(cmd_duohmm),
                                   "gzip {0} {1}".format(self.duohmm_haps,
                                                         self.duohmm_sample)],
                                   self.duohmm_haps + '.gz')
        self.settings['params'] = parse_command(parameters)

    def run_single(self, si_args, dm_args):

        qsub_parameters = ['-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.dirs['log']]
        # get checksum
        self.settings['checksum'] = utils.md5_for_file(self.checksum_file)
        yaml.dump(self.settings,
                  stream=open(self.param_f, 'w'))

        print sh.qsub('-N', self.run_id, qsub_parameters,
                      si_args, self.si_script)
        if self.duohmm_script is not None:
            print sh.qsub('-hold_jid', self.run_id,
                          '-N', 'duohmm' + self.run_id,
                          qsub_parameters + dm_args, self.duohmm_script)


# collection of methods for submission of MERLIN jobs with various defaults
# http://www.sph.umich.edu/csg/abecasis/merlin/
class Merlin(tool.Tool):

    @staticmethod
    def _get_version(executable):

        string = os.popen(executable).read()

        p = re.compile("MERLIN (.+) -")  # parentheses for capture groups
        m = p.search(str(string))
        if m:
            return m.group(1)
        else:
            print(string)
            raise Exception('Version not parsed.')

    def __init__(self, parameters, outdir, executable='merlin', version=None,
                 run_id=None):

        # edit parameters
        # no need to auto determine any parameters for merlin
        if version is None:
            version = self._get_version(executable)

        tool.Tool.__init__(self,
                           parameters=parameters,
                           executable=executable,
                           outdir=outdir,
                           name=Merlin.__name__,
                           version=version,
                           run_id=run_id)

    def parse_output(self):
        # optionally takes output dir
        pass
