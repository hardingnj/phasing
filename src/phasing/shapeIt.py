__author__ = 'Nicholas Harding'

import re
import os
import tool

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
        pass
        # (optionally can be passed a output directory)
        # returns a genotype matrix in numpy format, matching anhima specs
        # and a list with the sample name of each column

