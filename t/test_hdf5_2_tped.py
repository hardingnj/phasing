import unittest
import phasing.phasing
import phasing.hdf5_2_tped
import collections
import numpy as np
compare_list = lambda x, y: collections.Counter(x) == collections.Counter(y)

class convert_genotypes_to_strings(unittest.TestCase):
    def test_allref(self):
        genotypes_at_pos = np.zeros(20).reshape(10,2)
        self.assertTrue(compare_list(phasing.hdf5_2_tped.convert_gts_to_strings(genotypes_at_pos, 'A', 'G'), 
                        ["A A"]*10)) 

    def test_allalt(self):
        genotypes_at_pos = np.ones(20).reshape(10,2)
        self.assertTrue(compare_list(phasing.hdf5_2_tped.convert_gts_to_strings(genotypes_at_pos, 'A', 'G'), 
                        ["G G"]*10)) 
    
    def test_each(self):
        genotypes_at_pos = np.array([0, 0, 1, 1, 0, 1, -1, -1]).reshape(-1 ,2)
        self.assertTrue(compare_list(phasing.hdf5_2_tped.convert_gts_to_strings(genotypes_at_pos, 'A', 'G'), 
                        ["A A", "G G", "A G", "0 0"])) 

class return_correct_row(unittest.TestCase):
    def test_each(self):
        genotypes_at_pos = np.array([0, 0, 1, 1, 0, 1, -1, -1]).reshape(-1 ,2)
        self.assertEqual(phasing.hdf5_2_tped.get_tped_row(genotypes_at_pos, 'A', 'G', 1000001, '3L'),
                         "\t".join(['3L', 'snp1000001', '0', '1000001', "A A", "G G", "A G", "0 0"])) 

#def get_tped_row(gt_data, reference, alternate, position, contig):

