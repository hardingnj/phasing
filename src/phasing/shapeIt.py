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


class ShapeIt(tool.Tool):

    def _get_version(self):

        shapeit_v = os.popen(self.executable + ' --version').read()

        p = re.compile("Version : (.+)\n")  # parentheses for capture groups
        m = p.search(str(shapeit_v))
        if m:
            return m.group(1)
        else:
            print(shapeit_v)
            raise Exception('Version not parsed.')

    def __init__(self, parameters, outdir, executable='shapeit'):

        tool.Tool.__init__(self,
                           executable=executable,
                           outdir=outdir,
                           name='shapeIt',
                           version=self._get_version)

        # automatically work out output filenames.
        self.haplotypes_f = os.path.join(self.outdir, 'haplotypes')
        self.phased_f = os.path.join(self.outdir, 'phased')
        assert '--output-max' not in parameters

        # now we can call the parse method to init command string and dict
        self.parse_command(parameters + ['--output-max', self.haplotypes_f,
                                         self.phased_f])

    def parse_output(self):
        # (optionally can be passed a output directory) or at least
        # reconstructed
        # returns a genotype matrix in numpy format, matching anhima specs
        # and a list with the sample name of each column
        haplotype_data = pd.read_csv(self.haplotypes_f, sep=" ", header=None)

        haplotype_data = haplotype_data.rename(columns={
            0: 'contig', 1: 'id', 2: 'pos', 3: 'A1', 4: 'A2'})

        sample_data = pd.read_csv(self.phased_f, sep=" ", header=0,
                                  skiprows=0)
        sample_data = sample_data.ix[1:]

        # get ped genotypes into anhima numpy formats.
        htype_pos = np.array(haplotype_data.pos.values)

        # store as 3D genotypes
        htype_gts = np.array(haplotype_data.ix[:, 5:]).reshape((
            haplotype_data.shape[0], -1, 2))

        assert htype_gts.shape[1] == sample_data.shape[0]

        return htype_gts, sample_data.columns, {'pos': htype_pos}