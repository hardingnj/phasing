__author__ = 'Nicholas Harding'

import unittest
import phasing
import tempfile


class TestMerlin(unittest.TestCase):

    def test_initialize(self):

        tmp = tempfile.mkdtemp()
        params = ['--B', 'file', '--duoHMM']

        test_run = phasing.merlin.Merlin(params, tmp)
        self.assertEquals(test_run.version, '1.1.2')
        self.assertEquals(test_run.name, 'Merlin')