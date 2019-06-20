import unittest
from tcrdist import pairwise

class test_sequence_pair(unittest.TestCase):

    # Check that seqUence_pair enforces strings as inputs
    def test_sequence_pair_TypeError(self):
        with self.assertRaises(TypeError):
            pairwise.sequence_pair("CAGQASQGNLIF", 1)

    def test_leninput1_lt_5_implies_distance_is_none(self):
        self.assertTrue(pairwise.sequence_pair("CAGQASQGNLIF", "CAGA").hamming_distance == None)

    def test_leninput2_lt_5_implies_distance_is_none(self):
        self.assertTrue(pairwise.sequence_pair("CAGA", "CAGQASQGNLIF").hamming_distance == None)

    # THIS TEST WILL LIKELY BE MODIFIED TO BE LESS STRINGENT
    def test_lenalign_lt_min_leninput_implies_distance__is_none(self):
        self.assertTrue(pairwise.sequence_pair("CAGQASQGNLIF", "SSSSSSRR").hamming_distance == None)

    def test_sequence_pair_hamming_distance0(self):
        self.assertTrue(pairwise.sequence_pair("CAGQASQGNLIF", "CAGQASQGNLIF").hamming_distance == 0)

    def test_sequence_pair_hamming_distance0_083(self):
        self.assertTrue(round(pairwise.sequence_pair("CAGQASQGNLIF", "CAGQRSQGNLIF").hamming_distance, 3) == 0.083)


if __name__ == '__main__':
    unittest.main()
