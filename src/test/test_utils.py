__author__ = 'Nicholas Harding'


import unittest
import phasing
#import tempfile
#import random
import numpy as np


class TestMerlin(unittest.TestCase):

    def test_random(self):

        # create array of 1s and 2s.
        random_array = np.random.choice((1, 2), 10000).reshape((-1, 2))
        res_random = phasing.utils.calculate_switch_error(random_array)

        self.assertItemsEqual(res_random.shape, (2,))
        self.assertTrue(np.all(res_random > 0))

    def test_no_se(self):

        # create array of all 1s
        test_array = np.hstack([np.ones((10000, 2)), np.ones((10000, 2)) + 1])
        self.assertItemsEqual(test_array.shape, (10000, 4))

        res_test = phasing.utils.calculate_switch_error(test_array)
        self.assertItemsEqual(res_test, (0, 0, 0, 0))
