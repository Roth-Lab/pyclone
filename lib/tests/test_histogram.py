'''
Created on 2012-05-15

@author: andrew
'''
import unittest

import numpy as np
import random

from pyclone.utils import histogram

class Test(unittest.TestCase):
    def test_left_endpoint(self):
        values = [-1, ]
        
        counts = histogram(values, num_bins=2, min_value= -1, max_value=1)
        
        self.assertEqual(counts, [1, 0])

    def test_right_endpoint(self):
        values = [1, ]
        
        counts = histogram(values, num_bins=2, min_value= -1, max_value=1)
        
        self.assertEqual(counts, [0, 1])

    def test_counts_vs_numpy(self):
        bins = 100
        
        for left, right in [(-2, -1), (-0.5, 0.5), (1, 2), (100, 1000)]:    
            values = self._generate_random_values(left, right, 1000)
            
            counts = histogram(values, num_bins=bins, min_value=left, max_value=right)
            np_counts = list(np.histogram(values, bins=bins, range=(left, right))[0])
            
            self.assertEqual(np_counts, counts)
    
    def test_normalise_vs_numpy(self):
        bins = 100
        
        for left, right in [(-2, -1), (-0.5, 0.5), (1, 2), (100, 1000)]:    
            values = self._generate_random_values(left, right, 1000)
            
            counts = histogram(values, num_bins=bins, min_value=left, max_value=right, normalise=True)
            np_counts = list(np.histogram(values, bins=bins, range=(left, right), density=True)[0])
            
            for x, y in zip(counts, np_counts):
                self.assertAlmostEqual(x, y, 3)
    
    def test_values_less_than_min(self):
        values = [-0.01, ]
        
        counts = histogram(values, num_bins=2, min_value= 0, max_value=1)
        
        self.assertEqual(counts, [0, 0])
    
    def test_values_less_than_max(self):
        values = [1.01, ]
        
        counts = h`istogram(values, num_bins=2, min_value= 0, max_value=1)
        
        self.assertEqual(counts, [0, 0])
    
    def _generate_random_values(self, left, right, num):
        return [(right - left) * random.random() + left for _ in  range(num)]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
