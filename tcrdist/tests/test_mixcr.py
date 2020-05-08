import pytest
from tcrdist.repertoire import TCRrep
from tcrdist import mixcr
import numpy as np
import pandas as pd
import os

# INTEGRATION TESTS
def test_mixcr_to_tcrdist_on_clones():
    test_clones = os.path.join('tcrdist','test_files_compact','SRR5130260.1.test.fastq.output.clns.txt')
    df = mixcr.mixcr_to_tcrdist2(chain = "delta", organism = "human", clones_fn = test_clones)

    assert isinstance(df, pd.DataFrame)
    df1 = mixcr.remove_entries_with_invalid_vgene(df, chain = "delta", organism = "human")
    assert isinstance(df, pd.DataFrame)
    df1['subject'] = 'SRR5130260.1'
    
    tr = TCRrep(cell_df = df1, 
                organism = "human", 
                chains = ['delta'], 
                db_file='gammadelta_db.tsv')  
    print(tr.cell_df.shape[0])
    
    tr.infer_cdrs_from_v_gene(chain = 'delta', imgt_aligned=True)
    
    tr.index_cols = ['subject', "v_d_gene", 'd_d_gene', 'j_d_gene',
                    'cdr3_d_nucseq', 'cdr3_d_aa','cdr1_d_aa', 
                    'cdr2_d_aa', 'pmhc_d_aa']
    
    tr.deduplicate()
    assert isinstance(tr.clone_df, pd.DataFrame)

def test_mixcr_to_tcrdist_on_seqs():
    test_seqs = os.path.join('tcrdist','test_files_compact','SRR5130260.1.test.fastq.result.txt')
    df = mixcr.mixcr_to_tcrdist2(chain = "delta", organism = "human", seqs_fn = test_seqs, clones_fn = None)
    assert isinstance(df, pd.DataFrame)
    df1 = mixcr.remove_entries_with_invalid_vgene(df, chain = "delta", organism = "human")
    assert isinstance(df, pd.DataFrame)
    df1['subject'] = 'SRR5130260.1'
    df1['count'] = 1
    assert df1.shape[0] == 249
    tr = TCRrep(cell_df = df1, 
                organism = "human", 
                chains = ['delta'], 
                db_file='gammadelta_db.tsv')  
    print(tr.cell_df.columns)
    assert tr.cell_df.shape[0] == 249
    assert tr.cell_df.shape[0] == 249
    assert tr.cell_df.shape[1] == 7
    # ['v_d_gene', 'd_d_gene', 'j_d_gene', 'cdr3_d_nucseq', 'cdr3_d_aa','subject', 'count']
    tr.infer_cdrs_from_v_gene(chain = 'delta', imgt_aligned=True)
    assert tr.cell_df.shape[0] == 249
    assert tr.cell_df.shape[1] == 10
    # ['v_d_gene', 'd_d_gene', 'j_d_gene', 'cdr3_d_nucseq', 'cdr3_d_aa','subject', cdr1_d_aa', 'cdr2_d_aa', 'pmhc_d_aa']
    tr.index_cols = ['subject', "v_d_gene", 'd_d_gene', 'j_d_gene',
                    'cdr3_d_nucseq', 'cdr3_d_aa','cdr1_d_aa', 
                    'cdr2_d_aa', 'pmhc_d_aa']
    tr.deduplicate()
    assert tr.clone_df.shape[1] == 10


def test_mixcr_integration_with_correct_chain():
    test_clones_fn = os.path.join('tcrdist','test_files_compact','SRR5130260.1.test.fastq.output.clns.txt')

    df = mixcr.mixcr_to_tcrdist2(chain = "delta", organism = "human", seqs_fn = None, clones_fn = test_clones_fn )
    assert isinstance(df, pd.DataFrame)
    df1 = mixcr.remove_entries_with_invalid_vgene(df, chain = "delta", organism = "human")
    assert isinstance(df, pd.DataFrame)
    assert df1.shape[0] == 89

def test_mixcr_integration_with_wrong_chain():
    test_clones_fn = os.path.join('tcrdist','test_files_compact','SRR5130260.1.test.fastq.output.clns.txt')
    df = mixcr.mixcr_to_tcrdist2(chain = "gamma", organism = "human", seqs_fn = None, clones_fn = test_clones_fn )
    df2 = mixcr.remove_entries_with_invalid_vgene(df, chain = "gamma", organism = "human")
    assert df2.shape[0] == 0



# def test_mixcr_to_tcrdist_on_clones():
#     clones_seqs = os.path.join('tcrdist','test_files_compact','SRR5130260.1.test.fastq.result.txt')
#     df = mixcr.mixcr_to_tcrdist2(chain = "delta", organism = "human", seqs_fn = test_seqs, clones_fn = None)
#     assert isinstance(df, pd.DataFrame)
#     df1 = mixcr.remove_entries_with_invalid_vgene(df, chain = "delta", organism = "human")
#     assert isinstance(df, pd.DataFrame)
#     df1['subject'] = 'SRR5130260.1'
#     df1['count'] = 1
#     assert df1.shape[0] == 249
#     tr = TCRrep(cell_df = df1, 
#                 organism = "human", 
#                 chains = ['delta'], 
#                 db_file='gammadelta_db.tsv')  
#     print(tr.cell_df.columns)
#     assert tr.cell_df.shape[0] == 249
#     assert tr.cell_df.shape[0] == 249
#     assert tr.cell_df.shape[1] == 7
#     # ['v_d_gene', 'd_d_gene', 'j_d_gene', 'cdr3_d_nucseq', 'cdr3_d_aa','subject', 'count']
#     tr.infer_cdrs_from_v_gene(chain = 'delta', imgt_aligned=True)
#     assert tr.cell_df.shape[0] == 249
#     assert tr.cell_df.shape[1] == 10
#     # ['v_d_gene', 'd_d_gene', 'j_d_gene', 'cdr3_d_nucseq', 'cdr3_d_aa','subject', cdr1_d_aa', 'cdr2_d_aa', 'pmhc_d_aa']
#     tr.index_cols = ['subject', "v_d_gene", 'd_d_gene', 'j_d_gene',
#                     'cdr3_d_nucseq', 'cdr3_d_aa','cdr1_d_aa', 
#                     'cdr2_d_aa', 'pmhc_d_aa']
#     tr.deduplicate()
#     assert tr.clone_df.shape[1] == 10

    
# UNIT TESTS
def test_change_TRAVDV_to_TRAVdashDV_1():
    assert mixcr._change_TRAVDV_to_TRAVdashDV(s = 'TRAV29DV5*01') == 'TRAV29/DV5*01'
    assert mixcr._change_TRAVDV_to_TRAVdashDV(s = 'TRAV36DV7*01') == 'TRAV36/DV7*01'
    assert mixcr._change_TRAVDV_to_TRAVdashDV(s = 'TRAV36DV7*02') == 'TRAV36/DV7*02'
    
def test_change_TRAVDV_to_TRAVdashDV_2():
    assert mixcr._change_TRAVDV_to_TRAVdashDV(s='TRAV38-2DV8*01') == 'TRAV38-2/DV8*01'

def test_change_TRAVDV_to_TRAVdashDV_3():
    assert mixcr._change_TRAVDV_to_TRAVdashDV(s='TRAV38-1*01') == 'TRAV38-1*01'

def test_change_TRAVDV_to_TRAVdashDV_4():
    "NaN case"
    assert mixcr._change_TRAVDV_to_TRAVdashDV(s=np.NaN) is np.NaN 

def test_allele_00_to_01_1():
    assert mixcr._allele_00_to_01('TRDV3*00') == 'TRDV3*01'
    assert mixcr._allele_00_to_01('TRDV3*01') == 'TRDV3*01'
    assert mixcr._allele_00_to_01('TRDV3*02') == 'TRDV3*02'

def test_allele_00_to_01_2():
    assert mixcr._allele_00_to_01('TRDD3*00') == 'TRDD3*01'
    assert mixcr._allele_00_to_01('TRDD3*01') == 'TRDD3*01'
    assert mixcr._allele_00_to_01('TRDD3*02') == 'TRDD3*02'

def test_allele_00_to_01_3():
    assert mixcr._allele_00_to_01(np.NaN) is np.NaN

def test_take_top_mixcr_gene_hit():
    assert mixcr._take_top_mixcr_gene_hit('TRDD3*00(45),TRDD2*00(40)') == 'TRDD3*00'
    assert mixcr._take_top_mixcr_gene_hit('TRDD3*00(45)') == 'TRDD3*00'
    assert isinstance(mixcr._take_top_mixcr_gene_hit(np.NaN),float)
    assert mixcr._take_top_mixcr_gene_hit(np.NaN) is np.NaN

def test_validate_gene_names_1():
    df = pd.DataFrame({'v_d_gene':['TRDV3*01','TRDV1*01', 'TRAV29/DV5*01', 
                                  'TRAV38-2/DV8*01', "TRBV1*01"]})
    r = mixcr._validate_gene_names(  series = df['v_d_gene'], 
                                    chain = 'delta', 
                                    organism = 'human')
    assert np.all(r == pd.Series([True,True,True,True,False]))               

def test_validate_gene_names_2():
    df = pd.DataFrame({'v_d_gene':['TRDV3*01','TRDV1*01', 'TRAV29DV5*01', 
                                  'TRAV38-2DV8*01', "TRBV1*01"]})
    r = mixcr._validate_gene_names(  series = df['v_d_gene'], 
                                    chain = 'delta', 
                                    organism = 'human')
    assert np.all(r == pd.Series([True,True,False,False,False]))               