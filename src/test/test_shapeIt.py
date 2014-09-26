__author__ = 'Nicholas Harding'

import unittest
import phasing.algorithms
import tempfile


class TestShapeIt(unittest.TestCase):

    def test_initialize(self):

        tmp = tempfile.mkdtemp()
        params = ['--B', 'file', '--duohmm']

        test_run = phasing.algorithms.ShapeIt(params, tmp)

        self.assertIsInstance(test_run.tool_dict, dict)
        self.assertEquals(test_run.tool_dict['command']['--B'], 'file')
        self.assertTrue(test_run.tool_dict['command']['--duohmm'])

        self.assertEquals(test_run.tool_dict['command']['--output-max'],
                          test_run.haplotypes_f + ';' + test_run.phased_f)

    def test_numbers_ok(self):
        tmp = tempfile.mkdtemp()
        params = ['--thread', 0.05, '--setting', 1]

        t = phasing.algorithms.ShapeIt(params, tmp)
        self.assertEquals(t.tool_dict['command']['--thread'], '0.05')