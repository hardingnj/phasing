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

    def __init__(self, executable=None, outdir=None, name=None,
                 version=None, run_id=None, outfile=None):

        self.executable = executable
        self.command_string = ''
        self.command_dict = {}
        # nice to check that version is a string?
        if isinstance(version, basestring):
            self.version = version
        else:
            self.version = version()
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

    @classmethod
    def from_directory(cls, directory):

        # basename is run id
        directory = os.path.realpath(directory)
        assert os.path.isdir(directory)
        base, run_id = os.path.split(directory)

        param_yaml = os.path.join(directory, run_id + '_parameters.yaml')
        from_yaml = yaml.load(stream=open(param_yaml, 'r'))

        return cls(executable=from_yaml['executable'],
                   outdir=from_yaml['base_dir'],
                   name=from_yaml['name'],
                   version=from_yaml['version'],
                   run_id=run_id)

    def __str__(self):
        return "\n".join(['Outdir: ' + self.outdir,
                          'Executable: ' + self.executable,
                          yaml.dump(self.command_dict)])

    # GENERIC
    def parse_command(self, parameters):

        # create command
        command_dict = {'name': self.name,  'version': self.version,
                        'run_id': self.run_id, 'base_dir': self.basedir,
                        'executable': self.executable}
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

        self.command_string = ' '.join([self.executable] + [str(x)
                                                            for x in
                                                            parameters])
        self.command_dict = command_dict

    # GENERIC
    def dump_parameters(self):
        yaml.dump(self.command_dict,
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
