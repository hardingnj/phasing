__author__ = 'Nicholas Harding'

import os
import re
import yaml
import utils
import shutil
import sh
import uuid


# TO DO
# Add reconstruct run function, so that run is recreated from filepath
# Needs to detect whether shapeIt or Merlin etc. Class name can be added to
# parameters file
# checkout: http://stackoverflow.com/questions/141545/overloading-init-in-python


class Tool():

    def __init__(self, parameters=None, executable=None, outdir=None, name=None,
                 version=None, run_id=None, manipulate_parameters=None,
                 outfile=None):

        self.executable = executable
        self.version = version
        self.name = name
        self.outfile = outfile

        if run_id is None:
            self.run_id = name + '_' + str(uuid.uuid4().get_hex().upper()[0:8])
        else:
            self.run_id = run_id
        self.outdir = os.path.join(outdir, self.name, self.version, self.run_id)
        self.basedir = outdir

        # set log name and script name
        self.log_f = os.path.join(self.outdir, self.run_id + '.log')
        self.script_f = os.path.join(self.outdir, self.run_id + '.sh')
        self.param_f = os.path.join(self.outdir, self.run_id +
                                    '_parameters.yaml')

        # manipulate parameters happens here (if defined)
        if manipulate_parameters is not None:
            parameters = manipulate_parameters(parameters)
        self.tool_dict, self.command_string = self.parse_command(parameters)

    def __str__(self):
        return yaml.dump(self.tool_dict)

    # Method to create command and put all data in a yaml dict
    def parse_command(self, parameters):

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

        command_string = ' '.join([self.executable] + parameters)

        tool_dict = {'name': self.name,  'version': self.version,
                     'run_id': self.run_id, 'base_dir': self.basedir,
                     'executable': self.executable, 'command': command_dict,
                     'parameters': parameters}

        return tool_dict, command_string

    # GENERIC
    def dump_parameters(self):
        yaml.dump(self.tool_dict,
                  stream=open(self.param_f, 'w'))

    def delete(self):
        shutil.rmtree(self.outdir)

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
                               commands=['cd ' + self.outdir,
                                         self.command_string],
                               outfile=self.outfile)

        print sh.qsub(qsub_parameters + list(args) + [self.script_f])
