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

    def test_numbers_ok(self):
        tmp = tempfile.mkdtemp()
        params = ['--thread', 0.05, '--setting', 1]

        t = phasing.shapeIt.ShapeIt(params, tmp)
        self.assertEquals(t.command_dict['--thread'], '0.05')