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


class DuoHMM(tool.Tool):

    @staticmethod
    def _get_version(executable):

        duohmm_s = os.popen(executable).read()

        p = re.compile("Version: (.+)\n")  # parentheses for capture groups
        m = p.search(str(duohmm_s))
        if m:
            return m.group(1)
        else:
            print(duohmm_s)
            raise Exception('Version not parsed.')

    # This is a pretty stripped down class, just an add on to ShapeIt.
    def __init__(self, parameters, outdir, executable='duohmm', version=None,
                 run_id=None):
        if version is None:
            version = self._get_version(executable)

        tool.Tool.__init__(self,
                           parameters=parameters,
                           executable=executable,
                           outdir=outdir,
                           name=DuoHMM.__name__,
                           version=version,
                           run_id=run_id,
                           manipulate_parameters=self._manipulate_parameters)

    def _manipulate_parameters(self, parameters):

        #duohmm -H duo -O duohmm_haplotypes -G possible_GEs -R recomb_map
        self.inputroot = os.path.join(self.basedir, self.run_id)

        self.outfile = os.path.join(self.outdir, self.run_id)
        self.haplotypes_f = self.outfile + '.haps.gz'
        self.phased_f = self.outfile + '.sample'

        self.genotype_errors = os.path.join(self.outdir, 'genotype_errors.txt')
        self.recombination_map = os.path.join(self.outdir, 'recombination.txt')

        # looks for expected files, if not found renames them.
        # (For backwards compatability)
        if not os.path.isfile(os.path.join(self.basedir,
                                           self.run_id + '.haps')):
            try:
                os.rename(os.path.join(self.basedir, 'haplotypes'),
                          os.path.join(self.basedir, self.run_id + '.haps'))
                os.rename(os.path.join(self.basedir, 'phased'),
                          os.path.join(self.basedir, self.run_id + '.sample'))
            except OSError:
                pass

        for f in ('-H', '-O', '-G', '-R'):
            assert f not in parameters

        return parameters + ['-H', self.inputroot,
                             '-O', self.outfile,
                             '-G', self.genotype_errors,
                             '-R', self.recombination_map]


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
                 version=None):

        self.executable = executable
        self.ligate_bin = ligate_exec
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
        self.dirs = {d: os.path.join(self.outdir, d) for d in ('log', 'script')}
        [sh.mkdir(d, '-p') for d in self.dirs.values() if not os.path.isdir(d)]

        self.haplotypes_f = os.path.join(self.outdir, self.run_id +
                                         '_final.haps.gz')
        self.phased_f = os.path.join(self.outdir, self.run_id +
                                     '_final.samples.gz')
        self.param_f = os.path.join(self.outdir, self.run_id + '.yaml')
        self.settings = {'run_id': self.run_id, 'version': self.version,
                         'executable': self.executable, 'outdir': self.outdir,
                         'basedir': outdir, 'haps': self.haplotypes_f,
                         'samples': self.phased_f}

    def setup_region_jobs(self, parameters, regions=None, vcf_file=None,
                          pirs=None, duohmm=False):

        # basically, set up a shapeIT job for each region. This may be several
        # lines long.
        region_dir = os.path.join(self.outdir, 'region_outputs')
        os.mkdir(region_dir)

        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file

        hap_files = list()
        for i, region in enumerate(regions):
            start, stop = [str(x) for x in region]
            haps = os.path.join(self.outdir, region_dir, str(i) + '_' +
                                self.run_id + '.haps.gz')
            hap_files.append(haps)
            samples = os.path.join(self.outdir, region_dir, str(i) + '_' +
                                   self.run_id + '.sample.gz')

            cmd_shape_it = " ".join([self.executable] + parameters +
                                    ['--output-from', start, '--output-to',
                                     stop, '--output-max', haps, samples])

            script_name = os.path.join(self.dirs['script'], str(i) +
                                       '_shapeIt.sh')
            utils.create_sh_script(
                filename=script_name,
                commands=['cd ' + region_dir, cmd_shape_it],
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
        self.settings['params'] = parse_command(parameters)

    def parse_output(self):
        return ShapeIt.process_shapeit_output(self.haplotypes_f, self.phased_f)

    def run_region(self, si_args, li_args):

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

    @staticmethod
    def process_shapeit_output(haplotypes, samples):
        # returns a genotype matrix in numpy format, matching anhima specs
        # and a list with the sample name of each column
        genotype_cache = haplotypes + '_genotypes.npy'
        pos_cache = haplotypes + '_positions.npy'

        # load haplotype data:
        if os.path.isfile(genotype_cache):
            htype_gts = np.load(genotype_cache)
            htype_pos = np.load(pos_cache)
        else:
            haplotype_data = pd.read_csv(haplotypes, sep=" ", header=None)

            haplotype_data = haplotype_data.rename(columns={
                0: 'contig', 1: 'id', 2: 'pos', 3: 'A1', 4: 'A2'})
            # store as 3D genotypes
            htype_gts = np.array(haplotype_data.ix[:, 5:]).reshape((
                haplotype_data.shape[0], -1, 2))
            htype_pos = np.array(haplotype_data.pos.values)
            np.save(genotype_cache, htype_gts)
            np.save(pos_cache, htype_pos)

        sample_data = pd.read_csv(samples, sep=" ", header=0, skiprows=0)
        sample_data = sample_data.ix[1:]

        assert htype_gts.shape[1] == sample_data.shape[0]

        return htype_gts, sample_data.ID_2.tolist(), {'pos': htype_pos}


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
