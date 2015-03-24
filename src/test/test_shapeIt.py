#! /usr/bin/python

__author__ = 'Nicholas Harding'

import unittest
import phasing.algorithms
import tempfile
import tempfile
import os


class TestShapeIt(unittest.TestCase):

    def setUp(self):
        tmp = tempfile.NamedTemporaryFile(delete=False)
        with open(tmp.name + '.bed', 'wb') as fout:
            fout.write(os.urandom(1024))
        self.tempfile = tmp.name

    def test_initialize(self):

        tmp = tempfile.mkdtemp()
        params = ['-B', self.tempfile, '--duohmm']

        test_run = phasing.algorithms.ShapeIt(outdir=tmp)
        test_run.setup_region_jobs(params, duohmm=False,
                                   regions=[(1, 10), (11, 20)], pirs=None,
                                   vcf_file=self.tempfile)
        self.assertIsInstance(test_run.settings, dict)
        self.assertEquals(test_run.settings['params']['-B'], self.tempfile)
        self.assertTrue(test_run.settings['params']['--duohmm'])

    def test_numbers_ok(self):
        params = ['--thread', 0.05, '--setting', 1, '-B', self.tempfile]

        t = phasing.algorithms.ShapeIt(outdir='/tmp')

        t.setup_region_jobs(params, duohmm=False,
                            regions=[(1, 10), (11, 20)], pirs=None,
                            vcf_file=self.tempfile)
        self.assertEquals(t.settings['params']['--thread'], '0.05')

    def test_graphs(self):

        tmp = tempfile.mkdtemp()
        params = ['-B', self.tempfile, '--duohmm']

        test_run = phasing.algorithms.ShapeIt(outdir=tmp)
        test_run.setup_region_jobs(params, duohmm=False,
                                   regions=[(1, 10), (11, 20)], pirs=None,
                                   graphs=True, vcf_file=self.tempfile)

        self.assertIsInstance(test_run.settings, dict)
        self.assertEquals(test_run.settings['params']['-B'], self.tempfile)
        self.assertTrue(test_run.settings['params']['--duohmm'])

    def test_graphs_errors(self):
        # if PIRs specified should error
        tmp = tempfile.mkdtemp()
        params = ['-B', self.tempfile, '--duohmm']
        test_run = phasing.algorithms.ShapeIt(outdir=tmp)

        self.assertRaises(ValueError, test_run.setup_region_jobs, params,
                          duohmm=False, regions=[(1, 10), (11, 20)],
                          pirs='file', graphs=True, vcf_file=self.tempfile)