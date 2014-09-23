__author__ = 'Nicholas Harding'

import os
import re
import pprint
import  utils
import sh
import uuid


class Tool():

    def __init__(self, executable=None, outdir=None, name=None,
                 version=None, outfile=None):

        self.executable = executable
        self.command_string = ''
        self.command_dict = {}
        # nice to check that version is a string? 
        self.version = version()
        self.name = name
        self.outfile = outfile

        self.run_id = name + str(uuid.uuid4().get_hex().upper()[0:8])
        self.outdir = os.path.join(outdir, self.name, self.version, self.run_id)

        # set log name and script name
        self.log_f = os.path.join(self.outdir, self.run_id + '.log')
        self.script_f = os.path.join(self.outdir, self.run_id + '.sh')
        self.param_f = os.path.join(self.outdir, self.run_id +
                                    '_parameters.txt')

    # GENERIC
    def parse_command(self, parameters):

        # create command
        command_dict = {'version': self.version, 'run_id': self.run_id}  # later
        #  we can populate with defaults / allowed

        cl = re.compile('^-')
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

    # GENERIC
    def dump_parameters(self):
        pprint.pprint(self.command_dict,
                      stream=open(self.param_f, 'w'),
                      indent=2,
                      width=80,
                      depth=None)

    # GENERIC
    def run(self, *args):

        qsub_parameters = ['-N', self.run_id,
                           '-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.log_f]

        # kicks off qsub job with parameters
        os.makedirs(self.outdir)  # No check as should not exist already
        self.dump_parameters()

        utils.create_sh_script(filename=self.script_f,
                               commands=[self.command_string],
                               outfile=self.outfile)

        print sh.qsub(qsub_parameters + list(args) + [self.script_f])