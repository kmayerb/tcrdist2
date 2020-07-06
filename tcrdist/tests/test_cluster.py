import os
import pandas as pd
import pytest
import numpy as np

from tcrdist import cluster
from tcrsampler.sampler import TCRsampler
import collections

fn = os.path.join('tcrdist','tests','dash_pa_clone_df.csv')
test_clone_df = pd.read_csv(fn)

fn = os.path.join('tcrdist','tests','dash_pa_cluster_df.csv')
test_cluster_df = pd.read_csv(fn)

test_row = test_cluster_df.iloc[0,:]

ts = TCRsampler(default_background = "wirasinha_mouse_beta_t_1.tsv.sampler.tsv")

def test_that_inputs_exist():
    assert isinstance(test_clone_df, pd.DataFrame)
    assert isinstance(test_cluster_df, pd.DataFrame)

def test_row_get_cluster_attributes():
    test_row = test_cluster_df.iloc[0,:]
    print(test_row)
    cl_attrs = cluster._get_cluster_attributes(row = test_row, df= test_clone_df)

def test_get_cl_attributes_get_palmotif():
    test_row = test_cluster_df.iloc[0,:]
    print(test_row)
    cl_attrs = cluster._get_cluster_attributes(row = test_row, df= test_clone_df)
    assert isinstance(cl_attrs.neighbors, list)
    assert isinstance(cl_attrs.K_neighbors, int)
    assert isinstance(cl_attrs.ns_cdr3_seqs, list)
    assert isinstance(cl_attrs.ns_df, pd.DataFrame)
    assert isinstance(cl_attrs.ns_v_j_cdr3, pd.DataFrame)
    assert isinstance(cl_attrs.gene_usage, list)
    cl_motif = cluster._get_palmotifs(cluster_attributes = cl_attrs, sampler = ts, write = False)
    assert isinstance(cl_motif.raw_motif, tuple)
    assert isinstance(cl_motif.raw_motif[0], pd.DataFrame)
    assert isinstance(cl_motif.raw_motif[1], pd.DataFrame)
    assert isinstance(cl_motif.pal_motif, tuple)
    assert isinstance(cl_motif.pal_motif[0], pd.DataFrame)
    assert isinstance(cl_motif.pal_motif[1], pd.DataFrame)
    assert isinstance(cl_motif.bkgd_motif, tuple)
    assert isinstance(cl_motif.bkgd_motif[0], pd.DataFrame)
    assert isinstance(cl_motif.bkgd_motif[1], pd.DataFrame)

def test_get_cl_attributes_get_palmotif_and_write():
    test_row = test_cluster_df.iloc[0,:]
    print(test_row)
    cl_attrs = cluster._get_cluster_attributes(row = test_row, df= test_clone_df)
    cl_motif = cluster._get_palmotifs(cluster_attributes = cl_attrs, sampler = ts, write = True, dest = 'static')

def test_iterate_over_cluster_df_write_palmotifs():
    for i,row in test_cluster_df.head(10).iterrows():
        cl_attrs = cluster._get_cluster_attributes(row = row, df= test_clone_df)
        cl_motif = cluster._get_palmotifs(cluster_attributes = cl_attrs, sampler = ts, write = True, dest = 'static')

    contents = os.listdir('static3')
    for i in test_cluster_df.head(10)['cluster_id']:
        assert os.path.isfile(os.path.join('static3', f"{i}.pal_svg.svg"))
        assert os.path.isfile(os.path.join('static3', f"{i}.raw_svg.svg"))
        assert os.path.isfile(os.path.join('static3', f"{i}.bkgd_svg.svg"))


# def test_write_cluster_table():
#     test_row = test_cluster_df.iloc[0,:]
#     print(test_row)
#     cl_attrs = cluster._get_cluster_attributes(row = test_row, df= test_clone_df)
#     cluster._write_cluster_table(cl_attrs.ns_df, cluster_id= cl_attrs.cluster_id, cols = ['count'], chain = 'beta', dest='static')



def test__get_neighbors():
    r = cluster._get_neighbors(df = test_clone_df, 
                               ids = [16, 25, 29, 32, 61, 69, 94, 103, 108, 149, 163, 180, 183, 241, 265, 284, 323])
    assert isinstance(r, pd.DataFrame)
    assert r.shape[0] == len([16, 25, 29, 32, 61, 69, 94, 103, 108, 149, 163, 180, 183, 241, 265, 284, 323])
    assert np.all(r.index.to_list() == [16, 25, 29, 32, 61, 69, 94, 103, 108, 149, 163, 180, 183, 241, 265, 284, 323])

def test__get_neighbors_from_row():
    r = cluster._get_neighbors_from_row(row = test_row, df = test_clone_df)
    assert isinstance(r, pd.DataFrame)

def test__get_v_j_cdr3(chain = 'beta'):
    r = cluster._get_neighbors_from_row(row = test_row, df = test_clone_df)
    r2 = cluster._get_v_j_cdr3(df = r, chain = chain )
    r2b = cluster._get_v_j_cdr3(df = r, 
                            cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'], 
                            chain = 'beta')
    assert isinstance(r2, pd.DataFrame)
    assert np.all(r2 == r2b)







# # # INTEGRATION TESTS 
# # def test_get_cluster_attributes():
# #     c = cluster._get_cluster_attributes(row = test_row, df= test_clone_df)
# #     assert isinstance(c.ns_df, pd.DataFrame)
# #     assert isinstance(c.ns_v_j_cdr3, pd.DataFrame)
# #     assert isinstance(c.ns_cdr3_seqs, list)
# #     assert isinstance(c.gene_usage, list)
# #     assert isinstance(c.centroid, str)
# #     print(c.ns_cdr3_seqs)




#     # # setup a context
#     # import pandas as pd
#     # from tcrdist.repertoire import TCRrep
#     # from tcrdist import cluster
#     # import pwseqdist as pw
#     # df = pd.read_csv('dash.csv')
#     # df = df[df.epitope.isin(['PA'])]
#     # tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
#     # tr.tcrdist2(processes = 1,
#     #             metric = 'nw',
#     #             reduce = True,
#     #             dump = False,
#     #             save = False)

#     # from tcrsampler.sampler import TCRsampler
#     # t = TCRsampler(use_default=True)

#     # tr.cluster_index = tr.simple_cluster_index(pw_distances = tr.pw_beta,
#     #                                            method = 'complete',
#     #                                            criterion = "distance",
#     #                                            t = 75)

#     # assert len(tr.cluster_index) == tr.clone_df.shape[0]
#     # tr.cluster_df = tr.cluster_index_to_df(tr.cluster_index)


# # # INTEGRATION TESTS 
# # def test_get_cluster_attributes():
# #     c = cluster._get_cluster_attributes(row = test_row, df= test_clone_df)
# #     assert isinstance(c.ns_df, pd.DataFrame)
# #     assert isinstance(c.ns_v_j_cdr3, pd.DataFrame)
# #     assert isinstance(c.ns_cdr3_seqs, list)
# #     assert isinstance(c.gene_usage, list)
# #     assert isinstance(c.centroid, str)
# #     print(c.ns_cdr3_seqs)

# # # UNIT TESTS
# # def test__get_neighbors():
# #     r = cluster._get_neighbors(df = test_clone_df, 
# #                                ids = [534, 535, 541, 543, 1737, 1738, 2014, 2015, 2017])
# #     assert isinstance(r, pd.DataFrame)
# #     assert r.shape[0] == len([534, 535, 541, 543, 1737, 1738, 2014, 2015, 2017])
# #     assert np.all(r.index.to_list() == [534, 535, 541, 543, 1737, 1738, 2014, 2015, 2017])

# # def test__get_neighbors_from_row():
# #     r = cluster._get_neighbors_from_row(row = test_row, df = test_clone_df)
# #     assert isinstance(r, pd.DataFrame)

# # def test__get_v_j_cdr3(chain = 'beta'):
# #     r = cluster._get_neighbors_from_row(row = test_row, df = test_clone_df)
# #     r2 = cluster._get_v_j_cdr3(df = r, chain = chain )
# #     r2b = cluster._get_v_j_cdr3(df = r, 
# #                             cols = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'], 
# #                             chain = 'beta')
# #     assert isinstance(r2, pd.DataFrame)
# #     assert np.all(r2 == r2b)











# # import pytest
# # from tcrdist import cluster 
# # import pandas as pd
# # import numpy as np

# # test_row = pd.Series({'ct_columns': ('status', 'cmember'),
# #  'val_0': ('C', 'MEM+'),
# #  'ct_0': 9,
# #  'val_1': ('C', 'MEM-'),
# #  'ct_1': 2226,
# #  'val_2': ('H', 'MEM+'),
# #  'ct_2': 0,
# #  'val_3': ('H', 'MEM-'),
# #  'ct_3': 5367,
# #  'levels': [['C', 'H'], ['MEM+', 'MEM-']],
# #  'X+MEM+': 0,
# #  'X+MEM-': 5367,
# #  'X-MEM+': 9,
# #  'X-MEM-': 2226,
# #  'X_marg': 0.7059984214680347,
# #  'MEM_marg': 0.0011838989739542227,
# #  'X|MEM+': 0.0,
# #  'X|MEM-': 0.7068352429869617,
# #  'MEM|X+': 0.0,
# #  'MEM|X-': 0.004026845637583893,
# #  'index': 1738,
# #  'neighbors': [534, 535, 541, 543, 1737, 1738, 2014, 2015, 2017],
# #  'K_neighbors': 9,
# #  'R_radius': 50,
# #  'RR': 3.411051212938005,
# #  'OR': 100000,
# #  'pvalue': 1.622557654819239e-05,
# #  'FWERp': 0.12079941740129234,
# #  'FDRq': 0.000756729036315083})

# # test_clone_df = pd.DataFrame({'subject': {1: 1483,
#   534: 3211,
#   535: 3211,
#   541: 3211,
#   543: 3211,
#   1737: 4973,
#   1738: 4973,
#   2014: 6913,
#   2015: 6913,
#   2017: 6913},
#  'cell_type': {1: 'PBMC',
#   534: 'PBMC',
#   535: 'PBMC',
#   541: 'PBMC',
#   543: 'PBMC',
#   1737: 'PBMC',
#   1738: 'PBMC',
#   2014: 'PBMC',
#   2015: 'PBMC',
#   2017: 'PBMC'},
#  'v_b_gene': {1: 'TRBV10-3*01',
#   534: 'TRBV27*01',
#   535: 'TRBV27*01',
#   541: 'TRBV27*01',
#   543: 'TRBV27*01',
#   1737: 'TRBV27*01',
#   1738: 'TRBV27*01',
#   2014: 'TRBV27*01',
#   2015: 'TRBV27*01',
#   2017: 'TRBV27*01'},
#  'j_b_gene': {1: 'TRBJ2-7*01',
#   534: 'TRBJ2-4*01',
#   535: 'TRBJ2-4*01',
#   541: 'TRBJ2-4*01',
#   543: 'TRBJ2-4*01',
#   1737: 'TRBJ2-4*01',
#   1738: 'TRBJ2-4*01',
#   2014: 'TRBJ2-4*01',
#   2015: 'TRBJ2-4*01',
#   2017: 'TRBJ2-4*01'},
#  'cdr3_b_aa': {1: 'CAISESGSIHEQYF',
#   534: 'CASSALTGPAIAKNIQYF',
#   535: 'CASSKLTGEAVAKNIQYF',
#   541: 'CASSPLTGEALAKNIQYF',
#   543: 'CASSPLTGSGAAKNIQYF',
#   1737: 'CASSPLVGGPEAKNIQYF',
#   1738: 'CASSPLVGTALAKNIQYF',
#   2014: 'CASSPLTGALLAKNIQYF',
#   2015: 'CASSPLVGAPEAKNIQYF',
#   2017: 'CASSPSTGAAEAKNIQYF'},
#  'cdr2_b_aa': {1: 'SYG....VKD',
#   534: 'SMN....VEV',
#   535: 'SMN....VEV',
#   541: 'SMN....VEV',
#   543: 'SMN....VEV',
#   1737: 'SMN....VEV',
#   1738: 'SMN....VEV',
#   2014: 'SMN....VEV',
#   2015: 'SMN....VEV',
#   2017: 'SMN....VEV'},
#  'cdr1_b_aa': {1: 'ENH.......RY',
#   534: 'MNH.......EY',
#   535: 'MNH.......EY',
#   541: 'MNH.......EY',
#   543: 'MNH.......EY',
#   1737: 'MNH.......EY',
#   1738: 'MNH.......EY',
#   2014: 'MNH.......EY',
#   2015: 'MNH.......EY',
#   2017: 'MNH.......EY'},
#  'pmhc_b_aa': {1: 'S.KTED',
#   534: 'K.EKRN',
#   535: 'K.EKRN',
#   541: 'K.EKRN',
#   543: 'K.EKRN',
#   1737: 'K.EKRN',
#   1738: 'K.EKRN',
#   2014: 'K.EKRN',
#   2015: 'K.EKRN',
#   2017: 'K.EKRN'},
#  'cdr3_b_nuc': {1: 'CTGGAGTCCGCTACCAGCTCCCAGACATCTGTGTACTTCTGTGCCATCAGTGAGTCCGGCAGTATCCACGAGCAGTACTTCGGGCCG',
#   534: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTGCCCTAACGGGGCCAGCGATAGCCAAAAACATTCAGTACTTCGGCGCC',
#   535: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTAAATTGACAGGGGAGGCTGTTGCCAAAAACATTCAGTACTTCGGCGCC',
#   541: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCGCTGACAGGGGAGGCACTAGCCAAAAACATTCAGTACTTCGGCGCC',
#   543: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCCCTGACAGGGTCAGGTGCGGCCAAAAACATTCAGTACTTCGGCGCC',
#   1737: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCCCTCGTCGGGGGCCCTGAGGCCAAAAACATTCAGTACTTCGGCGCC',
#   1738: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCCCTTGTTGGGACAGCCCTAGCCAAAAACATTCAGTACTTCGGCGCC',
#   2014: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCGTTGACAGGGGCCTTACTAGCCAAAAACATTCAGTACTTCGGCGCC',
#   2015: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCCCTGGTGGGAGCCCCCGAGGCCAAAAACATTCAGTACTTCGGCGCC',
#   2017: 'AGCCCCAACCAGACCTCTCTGTACTTCTGTGCCAGCAGTCCTTCCACAGGGGCCGCCGAGGCCAAAAACATTCAGTACTTCGGCGCC'},
#  'epitope': {1: 'AFPFTIYSL,GYINVFAFPF,INVFAFPFTI,MGYINVFAF,NVFAFPFTI,NVFAFPFTIY,YINVFAFPF',
#   534: 'LSPRWYFYY,SPRWYFYYL',
#   535: 'LSPRWYFYY,SPRWYFYYL',
#   541: 'LSPRWYFYY,SPRWYFYYL',
#   543: 'LSPRWYFYY,SPRWYFYYL',
#   1737: 'LSPRWYFYY,SPRWYFYYL',
#   1738: 'LSPRWYFYY,SPRWYFYYL',
#   2014: 'LSPRWYFYY,SPRWYFYYL',
#   2015: 'LSPRWYFYY,SPRWYFYYL',
#   2017: 'LSPRWYFYY,SPRWYFYYL'},
#  'count': {1: 1,
#   534: 1,
#   535: 1,
#   541: 1,
#   543: 1,
#   1737: 1,
#   1738: 1,
#   2014: 1,
#   2015: 1,
#   2017: 1},
#  'clone_id': {1: 2,
#   534: 535,
#   535: 536,
#   541: 542,
#   543: 544,
#   1737: 1738,
#   1738: 1739,
#   2014: 2015,
#   2015: 2016,
#   2017: 2018},
#  'status': {1: 'H',
#   534: 'C',
#   535: 'C',
#   541: 'C',
#   543: 'C',
#   1737: 'C',
#   1738: 'C',
#   2014: 'C',
#   2015: 'C',
#   2017: 'C'}})






