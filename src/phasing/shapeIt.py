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


class ShapeIt:

    def _get_version(self):

        shapeit_v = os.popen(self.executable + ' --version').read()

        p = re.compile("Version : (.+)\n")  # parentheses for capture groups
        m = p.search(str(shapeit_v))
        if m:
            self.version = m.group(1)
        else:
            print(shapeit_v)
            raise Exception('Version not parsed.')

    def __init__(self, parameters, outdir, executable='shapeit'):
        self.executable = executable
        self._get_version()
        self.command_string = ''
        self.command_dict = {}

        self.run_id = 'shapeIt_' + str(uuid.uuid4().get_hex().upper()[0:8])
        self.outdir = os.path.join(outdir, 'shapeit', self.version, self.run_id)

        # automatically work out output filenames.
        self.haplotypes_f = os.path.join(self.outdir, 'haplotypes')
        self.phased_f = os.path.join(self.outdir, 'phased')

        assert '--output-max' not in parameters
        self.parse_command(parameters + ['--output-max', self.haplotypes_f,
                                         self.phased_f])

        self.log_f = os.path.join(self.outdir, self.run_id + '.log')
        self.script_f = os.path.join(self.outdir, self.run_id + '.sh')

    def parse_command(self, parameters):

        # create command
        command_dict = {'version': self.version, 'run_id': self.run_id}  # later
        #  we can populate with defaults / allowed

        cl = re.compile('^--')
        last_was_key = False
        key = None

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

        self.command_string = ' '.join([self.executable] + parameters)
        self.command_dict = command_dict

    def dump_parameters(self):
        pprint.pprint(self.command_dict,
                      stream=open(self.log_f, 'w'),
                      indent=2,
                      width=80,
                      depth=None)

    def run(self, *args):

        qsub_parameters = ['-N', self.run_id,
                           '-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.log_f]

        # kicks off qsub job with parameters
        os.makedirs(self.outdir)  # No check as should not exist already
        self.dump_parameters()

        phasing.create_sh_script(self.script_f, [self.command_string],
                                 self.haplotypes_f)

        sh.qsub(qsub_parameters + list(args) + [self.script_f])