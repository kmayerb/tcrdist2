import unittest
import tcrdist
from tcrdist import repertoire_db
import numpy as np



class test_repertoire_db(unittest.TestCase):

    def test_extract_cdrs_without_gaps(self):
        """ specific case """
        a = repertoire_db._extract_cdrs_without_gaps('GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..',
     '28-39;57-66;82-88;106-111')
        b = ['TSGFNG', 'VVLDGL', 'SRSNGY', 'CAVR']
        self.assertTrue(a==b)

    def test_extract_aligned_cdrs(self):
        """ specific case """
        a = repertoire_db._extract_aligned_cdrs('GQGVEQ.P.AKLMSVEGTFARVNCTYSTSG......FNGLSWYQQREGQAPVFLSYVVL....DGLKDS.....GHFSTFLSRSN.GYSYLLLTELQIKDSASYLCAVR..',
     '28-39;57-66;82-88;106-111')
        b = ['TSG......FNG', 'VVL....DGL', 'SRSN.GY', 'CAVR..']
        self.assertTrue(a==b)

    def test_repertoire_db_extract_mouse_cdrs(self):
        """ tests that repertoire_db.extract_aligned_cdrs() is producing
        exactly the same hundred mouse results as are present in the column 'cdrs'
        in alphabeta_db.tsv """
        d = repertoire_db.generate_dict_from_db("alphabeta_db.tsv")
        mouse_cdrs_in_file = [d['mouse'][k]["cdrs"] for k in d['mouse'].keys()]
        mouse_cdrs_in_file = [x.split(';') for x in mouse_cdrs_in_file[0:100]]
        aligned_protseqs = [d['mouse'][k]["aligned_protseq"] for k in d['mouse'].keys()]
        cdr_columns = [d['mouse'][k]["cdr_columns"] for k in d['mouse'].keys()]
        mouse_cdrs_infered = list(map(lambda a,p: repertoire_db._extract_aligned_cdrs(a,p),
                                      aligned_protseqs[:], cdr_columns[:]))
        self.assertTrue(mouse_cdrs_infered[0:100] ==  mouse_cdrs_in_file)
        #assert(all([x == y for x, y in zip(mouse_cdrs_infered[0:100], mouse_cdrs_in_file)]))

    def test_repertoire_db_extract_human_cdrs(self):
        """ tests that repertoire_db.extract_aligned_cdrs() is producing
        exactly the same hundred human results as are present in the column 'cdrs'
        in alphabeta_db.tsv """

        d = repertoire_db.generate_dict_from_db("alphabeta_db.tsv")
        human_cdrs_in_file = [d['human'][k]["cdrs"] for k in d['human'].keys()]
        human_cdrs_in_file = [x.split(';') for x in human_cdrs_in_file[0:100]]
        aligned_protseqs = [d['human'][k]["aligned_protseq"] for k in d['human'].keys()]
        cdr_columns = [d['human'][k]["cdr_columns"] for k in d['human'].keys()]
        human_cdrs_infered = list(map(lambda a,p: repertoire_db._extract_aligned_cdrs(a,p),
                                      aligned_protseqs[:], cdr_columns[:]))
        self.assertTrue(human_cdrs_infered[0:100] ==  human_cdrs_in_file)
        #assert(all([x == y for x, y in zip(human_cdrs_infered[0:100], human_cdrs_in_file)]))






if __name__ == '__main__':
    unittest.main()
