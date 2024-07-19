#There is currently an issue whereby transcripts associated with our control gene (Limk1) are not showing up
#in the output of intersects.py. We therefore need to test if the script is working properly.

#There are a few different ways bf.overlap could be returning data:
#1: full overlap- both start and end of binding sites within UTR.
#2: full overlap 2- neither start or end of binding site within UTR, but UTR is within the binding site.
#3: partial overlap: only start coord in the UTR.
#4: partial overlap: only end coord in the UTR.
#5: no overlap: neither start nor end of binding site within UTR- UTR not within binding site.
#Start with the above.

#TEST_INPUTS
#3UTR_test = input test 3UTR coordinates
#DDX6_test = input test DDX6 binding site coordinates (reflects two examples of each of the above 1-5)
#TEST_OUTPUTS
#looks like the bioframe.overlap() function returns all full and partial overlaps.
#all_full_partial_overlaps = first possible output (all full + partial overlaps)
#full_overlaps_only = second possible output (all full overlaps only)
#partial_overlaps_only = third possible output (all partial overlaps)

#Ultimately two things to do: write a function that wraps around the entirety of intersects.py to test that 
#it works, and test the individual functions that I've written.

"""Tests for get_additional_columns and get_GENCODE_ID"""
#run the below tests using: python -m unittest test.test_module in the terminal

#get_additional_columns looks up additional information related to genomic coordinates from the DDX6_coords 
#DataFrame based on the input row. If a match is found, it returns additional info from that row; otherwise, it returns None.


import unittest
import pandas as pd
import numpy as np

import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from intersects import get_additional_columns


class TestGetAdditionalColumns(unittest.TestCase):
    def setUp(self):
        # Create a sample DataFrame for testing
        self.row = pd.DataFrame({
            'chrom': ['chr1'],
            'start': [100],
            'end': [200]
        })

        # Create a mock DDX6_coords DataFrame (replace with your actual data)
        self.DDX6_coords = pd.DataFrame({
            'chrom': ['chr1'],
            'start': [100],
            'end': [200],
            'strand': ['+'],
            'confidence_score': [60.5],
            'sample_or_tissue_used': ['K562']
        })

    def test_matching_row(self):
        # Test when a matching row exists
        result = get_additional_columns(self.row)
        expected_result = {
            'strand': '+',
            'confidence_score': 60.5,
            'sample_or_tissue_used': 'K562'
        }
        np.testing.assert_equal(result, expected_result)

    def test_no_matching_row(self):
        # Test when no matching row exists
        non_matching_row = pd.DataFrame({
            'chrom': ['chr2'],
            'start': [300],
            'end': [400]
        })
        result = get_additional_columns(non_matching_row)
        self.assertIsNone(result)

if __name__ == '__main__':
    unittest.main()
