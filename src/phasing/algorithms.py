__author__ = 'Nicholas Harding'

import re
import os
import tool
import pandas as pd
import numpy as np
# collection of methods for submission of shapeit jobs with various defaults

# behaviour:

# constructor will work out version
# print out handy yaml file with settings etc
# the class represents a run and keeps track of all the files etc.

# main function returns a string with the shapeit command

# also have a method that submits a qsub. This is generic

# methods that create the dir structure etc. Probably generic too

# data/shapeit/version/


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
        self.haplotypes_f = self.outfile + '.haps'
        self.phased_f = self.outfile + '.sample'

        self.genotype_errors = os.path.join(self.outdir, 'genotype_errors.txt')
        self.recombination_map = os.path.join(self.outdir, 'recombination.txt')

        # looks for expected files, if not found renames them.
        # (For backwards compatability)
        if not os.path.isfile(os.path.join(self.basedir,
                                           self.run_id + '.haps')):
            os.rename(os.path.join(self.basedir, 'haplotypes'),
                      os.path.join(self.basedir, self.run_id + '.haps'))
            os.rename(os.path.join(self.basedir, 'phased'),
                      os.path.join(self.basedir, self.run_id + '.sample'))

        for f in ('-H', '-O', '-G', '-R'):
            assert f not in parameters

        return parameters + ['-H', self.inputroot,
                             '-O', self.outfile,
                             '-G', self.genotype_errors,
                             '-R', self.recombination_map]

    def parse_output(self):
        return ShapeIt.process_shapeit_output(self.haplotypes_f,
                                              self.phased_f)


class ShapeIt(tool.Tool):

    @staticmethod
    def _get_version(executable):

        shapeit_v = os.popen(executable + ' --version').read()

        p = re.compile("Version : (.+)\n")  # parentheses for capture groups
        m = p.search(str(shapeit_v))
        if m:
            return m.group(1)
        else:
            print(shapeit_v)
            raise Exception('Version not parsed.')

    def _manipulate_parameters(self, parameters):
        # function edits parameters before they are passed to Tool
        # handy for auto determining outputs etc
        # automatically work out output filenames.
        self.haplotypes_f = os.path.join(self.outdir, self.run_id + '.haps')
        self.phased_f = os.path.join(self.outdir, self.run_id + '.sample')
        assert '--output-max' not in parameters

        # now we can call the parse method to init command string and dict
        return parameters + ['--output-max', self.haplotypes_f, self.phased_f]

    def __init__(self, parameters, outdir, executable='shapeit',
                 version=None, run_id=None):

        if version is None:
            version = self._get_version(executable)

        tool.Tool.__init__(self,
                           parameters=parameters,
                           executable=executable,
                           outdir=outdir,
                           name=ShapeIt.__name__,
                           version=version,
                           run_id=run_id,
                           manipulate_parameters=self._manipulate_parameters)

    def parse_output(self):
        return ShapeIt.process_shapeit_output(self.haplotypes_f, self.phased_f)

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
