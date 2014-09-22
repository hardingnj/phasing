__author__ = 'Nicholas Harding'


import re
import sh
import os
import uuid
import pprint
import phasing

# collection of methods for submission of shapeit jobs with various defaults

# behaviour:

# constructor will work out version
# print out handy yaml file with settings etc
# the class represents a run and keeps track of all the files etc.

# main function returns a string with the shapeit command

# also have a method that submits a qsub. This is generic

# methods that create the dir structure etc. Probably generic too

# data/shapeit/version/

# have a get_version method


class Merlin:

    def _get_version(self):

        string = os.popen(self.executable).read()

        p = re.compile("MERLIN (.+) -")  # parentheses for capture groups
        m = p.search(str(string))
        if m:
            self.version = m.group(1)
        else:
            print(string)
            raise Exception('Version not parsed.')

    def __init__(self, parameters, outdir, executable='merlin'):
        self.executable = executable
        self._get_version()
        self.command_string = ''
        self.command_dict = {}

        self.run_id = 'merlin_' + str(uuid.uuid4().get_hex().upper()[0:8])
        self.outdir = os.path.join(outdir, 'merlin', self.version, self.run_id)

        # add output names here

        # create log and script file
        self.log_f = os.path.join(self.outdir, self.run_id + '.log')
        self.script_f = os.path.join(self.outdir, self.run_id + '.sh')

    def parse_output(self):
        # optionally takes output dir
        pass