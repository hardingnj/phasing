__author__ = 'Nicholas Harding'

import re
import os
import tool


# collection of methods for submission of MERLIN jobs with various defaults
# http://www.sph.umich.edu/csg/abecasis/merlin/
class Merlin(tool.Tool):

    def _get_version(self):

        string = os.popen(self.executable).read()

        p = re.compile("MERLIN (.+) -")  # parentheses for capture groups
        m = p.search(str(string))
        if m:
            return m.group(1)
        else:
            print(string)
            raise Exception('Version not parsed.')

    def __init__(self, parameters, outdir, executable='merlin'):

        tool.Tool.__init__(self,
                           executable=executable,
                           outdir=outdir,
                           name='Merlin',
                           version=self._get_version)

        # no need to auto determine any parameters for merlin
        self.parse_command(parameters)

    def parse_output(self):
        # optionally takes output dir
        pass