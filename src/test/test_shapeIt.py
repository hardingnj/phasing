__author__ = 'Nicholas Harding'

import unittest
import phasing.shapeIt
import tempfile


class TestShapeIt(unittest.TestCase):

    def test_initialize(self):

        tmp = tempfile.mkdtemp()
        params = ['--B', 'file', '--duoHMM']

        test_run = phasing.shapeIt.ShapeIt(params, tmp)

        self.assertIsInstance(test_run.command_dict, dict)
        self.assertEquals(test_run.command_dict['--B'], 'file')
        self.assertTrue(test_run.command_dict['--duoHMM'])

        self.assertEquals(test_run.command_dict['--output-max'],
                          test_run.haplotypes_f + ';' + test_run.phased_f)

        #test_run.run('-l', 'h_vmem=2G')