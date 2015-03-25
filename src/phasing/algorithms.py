__author__ = 'Nicholas Harding'

import re
import os
import tool
import sh
import uuid
import tempfile
import utils
import subprocess
import yaml
import sys
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
