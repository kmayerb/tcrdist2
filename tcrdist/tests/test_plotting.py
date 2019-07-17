import unittest
import tcrdist
import pandas as pd
import numpy as np

class test_repertoire(unittest.TestCase):

    def test_gene_pairing_plot(self):
        np.random.seed(110820)
        n = 50
        df = pd.DataFrame({'VA':np.random.choice(['TRAV14', 'TRAV12', 'TRAV3', 'TRAV23', 'TRAV11', 'TRAV6'], n),
                           'JA':np.random.choice(['TRAJ4', 'TRAJ2', 'TRAJ3','TRAJ5', 'TRAJ21', 'TRAJ13'], n),
                           'VB':np.random.choice(['TRBV14', 'TRBV12', 'TRBV3', 'TRBV23', 'TRBV11', 'TRBV6'], n),
                           'JB':np.random.choice(['TRBJ4', 'TRBJ2', 'TRBJ3','TRBJ5', 'TRBJ21', 'TRBJ13'], n)})
        df = df.assign(Count=1)
        df.loc[:10, 'Count'] = 10
        svg = tcrdist.plotting.plotPairings(df, ['JA', 'VA', 'VB', 'JB'], count_col='Count')
        self.assertTrue(len(svg) > 0)

if __name__ == '__main__':
    unittest.main()
