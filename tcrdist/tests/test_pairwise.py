import unittest
from tcrdist import pairwise

class test_SequencePair(unittest.TestCase):

    # Tests gaps are introduced at sequence ends
    def test_align___0(self):
        # initialize SequencePair
        sp = pairwise.SequencePair("", "")
        # call align function with specific input
        self.assertTrue(sp.align("AAA", "GAAAGGG", open_penalty = 3, extend_penalty = 3).traceback.query == '-AAA---')

    # Tests that the complete longer sequence is returned
    def test_align___1(self):
        self.assertTrue(pairwise.SequencePair("AACAGQASQGNLIF", "CAGQASQGNLIF").a1 == 'AACAGQASQGNLIF')

    # Tests that the shorter sequence is returned with gaps
    def test_align___2(self):
        self.assertTrue(pairwise.SequencePair("AACAGQASQGNLIF", "CAGQASQGNLIF").a2 == '--CAGQASQGNLIF')
    # Check that SequencePair enforces strings as inputs
    def test_SequencePair___returns_TypeError(self):
        with self.assertRaises(TypeError):
            pairwise.SequencePair("CAGQASQGNLIF", 1)

    def test_SequencePair___leninput1_lt_5_returns_distance_None(self):
        self.assertTrue(pairwise.SequencePair("CAGA", "CAGQASQGNLIF").hamming_distance == None)

    def test_SequencePair___leninput2_lt_5_returns_distance_None(self):
        self.assertTrue(pairwise.SequencePair("CAGQASQGNLIF", "CAGA").hamming_distance == None)

    def test_SequencePair___returns_hamming_distance0(self):
        self.assertTrue(pairwise.SequencePair("CAGQASQGNLIF", "CAGQASQGNLIF").hamming_distance == 0)

    def test_SequencePair___returns_hamming_distance0_083(self):
        self.assertTrue(pairwise.SequencePair("CAGQASQGNLIF", "CAGQRSQGNLIF").hamming_distance== 1.0)

    def test_SequencePair___returns_hamming_distance0_5(self):
        self.assertTrue(pairwise.SequencePair("AAAAAAAAAA", "AAAAA").hamming_distance == 5.0)


    def test_unpack_dd_to_kkv(self):
        dd = { 'A': {'A': 1, 'B': 2, 'C': 3},
               'B': {'B': 4, 'C': 5},
               'C': {'C': 6}
             }
        expected = {'key1': ['A', 'A', 'A', 'B', 'B', 'C'],
        'key2': ['A', 'B', 'C', 'B', 'C', 'C'],
        'value': [1, 2, 3, 4, 5, 6]}

        kkv = pairwise.unpack_dd_to_kkv(dd)
        self.assertTrue(kkv == expected)

    def test_unpack_pooled_dd_to_kkv(self):
        pooled_dd = [
        { 'A': {'A': 1, 'B': 2, 'C': 3},
          'B': {'B': 4, 'C': 5},
          'C': {'C': 6}},
        { 'D': {'D': 7, 'E': 8, 'F': 9},
          'E': {'E': 10, 'F': 11},
          'F': {'F': 12}}
        ]

        expected  = {
        'key1'  : ['A', 'A', 'A', 'B', 'B', 'C', 'D', 'D', 'D', 'E', 'E', 'F'],
        'key2'  : ['A', 'B', 'C', 'B', 'C', 'C', 'D', 'E', 'F', 'E', 'F', 'F'],
        'value' : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        }

        kkv = pairwise.unpack_pooled_dd_to_kkv(pooled_dd)
        self.assertTrue(kkv == expected)


    def test_apply_pairwise_distance___returns_hamming_values(self):
        sequences = ["CAGQASQGNLIF", "CAGQASQGNLIF", "CAGQASQGNLIFA", \
        "CAGQASQGNLIAA", "CAGQASQGNLIFAAA", "CAGQASQGNLIFAAAAA"]
        # define it explicitly here since this is test for hamming values
        def my_distance_wrapper(a, b):
            return(float(pairwise.SequencePair(a, b).hamming_distance) )
        unique_seqs = pairwise.select_unique_sequences(sequences)
        d = pairwise.apply_pairwise_distance(unique_seqs)
        kkv = pairwise.unpack_dd_to_kkv(dd = d)
        a = kkv['value'][0:5]
        b = [0.0, 2.0, 4.0, 2.0, 0.0]
        self.assertTrue(a==b)

        # kkv- MANUAL INSPECTION
        # 'CAGQASQGNLIAA'
        # 'CAGQASQGNLIAA'  #0

        # 'CAGQASQGNLIAA'
        # 'CAGQASQGNLIFAAA' # 2

        # 'CAGQASQGNLIAA',
        # 'CAGQASQGNLIFAAAAA', # 4

        # 'CAGQASQGNLIF'
        # 'CAGQASQGNLIAA', # 2

        # 'CAGQASQGNLIF',
        # 'CAGQASQGNLIF', # 0

        #self.assertTrue(map(lambda p: round(p, 3), a) == map(lambda p: round(p, 3), b) )

    def test_apply_pairwise_distance___returns_correct_key1(self):

        sequences = ["CAGQASQGNLIF", "CAGQASQGNLIF", "CAGQASQGNLIFA", \
        "CAGQASQGNLIAA", "CAGQASQGNLIFAAA", "CAGQASQGNLIFAAAAA"]

        unique_seqs = pairwise.select_unique_sequences(sequences)
        d = pairwise.apply_pairwise_distance(unique_seqs)
        kkv = pairwise.unpack_dd_to_kkv(dd = d)
        a = kkv['key1'][0:5]
        b =['CAGQASQGNLIAA', 'CAGQASQGNLIAA', 'CAGQASQGNLIAA', 'CAGQASQGNLIF', 'CAGQASQGNLIF']
        self.assertTrue(a == b)

    def test_apply_pairwise_distance___gets_correct_key2(self):
        sequences = ["CAGQASQGNLIF", "CAGQASQGNLIF", "CAGQASQGNLIFA", \
        "CAGQASQGNLIAA", "CAGQASQGNLIFAAA", "CAGQASQGNLIFAAAAA"]
        unique_seqs = pairwise.select_unique_sequences(sequences)
        d = pairwise.apply_pairwise_distance(unique_seqs)
        kkv = pairwise.unpack_dd_to_kkv(dd = d)
        a = kkv['key2'][0:5]
        b =['CAGQASQGNLIAA', 'CAGQASQGNLIFAAA', 'CAGQASQGNLIFAAAAA', 'CAGQASQGNLIAA', 'CAGQASQGNLIF']
        self.assertTrue(a == b)

    #def test_apply_pairwise_distance_multiprocessing(self):
    #    sequences = ["CAGQASQGNLIF","CAGQASQGNLIF","CAGQASQGNLIFA", \
    #    "CAGQASQGNLIAA","CAGQASQGNLIFAAA","CAGQASQGNLIFAAAAA", "CAGQASQGNLIFG"]

    #    unique_seqs = pairwise.select_unique_sequences(sequences)
    #    d = pairwise.apply_pairwise_distance_multiprocessing(sequences = unique_seqs,  processes = 1)
    #    kkv = pairwise.unpack_pooled_dd_to_kkv(pooled_dd = d)
        #a = kkv['key2'][0:5]
        #b =['CAGQASQGNLIAA','CAGQASQGNLIFAAA','CAGQASQGNLIFAAAAA','CAGQASQGNLIAA','CAGQASQGNLIF']
        #a = kkv['key2'][0:5]
        # =['CAGQASQGNLIAA','CAGQASQGNLIFAAA','CAGQASQGNLIFAAAAA','CAGQASQGNLIAA','CAGQASQGNLIF']
    #    self.assertTrue(kkv)


if __name__ == '__main__':
    unittest.main()
