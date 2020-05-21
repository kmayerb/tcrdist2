import pytest
import pandas as pd
from tcrdist.subset import TCRsubset
from tcrdist.tests.my_test_subset import dist_d_subset, clone_df_subset_d
from tcrdist.tests.my_test_subset import dist_a_subset, dist_b_subset, clone_df_subset 
import numpy as np
from tcrdist.cdr3_motif import TCRMotif
from tcrdist import paths 

assert isinstance(dist_a_subset, pd.DataFrame)
assert isinstance(dist_b_subset, pd.DataFrame)
assert isinstance(clone_df_subset, pd.DataFrame)

def test_ng_tcrs_AB():
    """
    Test new method for ng_tcrs generation
    """
    df = clone_df_subset.iloc[0:20, :].copy()
    db = dist_b_subset.iloc[0:20, 0:20]
    da = dist_a_subset.iloc[0:20, 0:20]
    ts = TCRsubset(clone_df = df,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["A","B"],
            dist_a = da,
            dist_b = db)
    tm = TCRMotif(  clones_df = ts.tcr_motif_clones_df(), 
                    organism = ts.organism, 
                    chains = ts.chains, 
                    epitopes = ts.epitopes,
                    db_file = "alphabeta_db.tsv")
    ng_tcrs = dict()
    # Default Behavior
    for chain in tm.chains: #['A','B']:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref

    assert isinstance(ng_tcrs['A'], dict) 
    assert set(ng_tcrs['A'].keys()) == set(['TRAV14D-3/DV8*02', 'TRAV14-1*02', 'TRAV14-2*01', 'TRAV14D-1*01', 'TRAV14D-2*01', 'TRAV14D-2*02', 'TRAV14-3*01', 'TRAV14-3*02', 'TRAV14D-3/DV8*01', 'TRAV6-7/DV9*02', 'TRAV6-6*01', 'TRAV6-7/DV9*01', 'TRAV6-7/DV9*04', 'TRAV6D-6*01', 'TRAV6-5*01', 'TRAV6D-7*01', 'TRAV13-1*01', 'TRAV13D-4*01', 'TRAV13-3*01', 'TRAV10*01', 'TRAV2*01', 'TRAV6-4*01', 'TRAV3D-3*02', 'TRAV3-1*01', 'TRAV3-4*01', 'TRAV4D-4*03', 'TRAV4D-3*01', 'TRAV4D-4*01', 'TRAV4-3*01', 'TRAV4-2*01', 'TRAV4-4/DV10*01', 'TRAV4D-3*03', 'TRAV5-4*01', 'TRAV5D-4*02', 'TRAV3-3*01', 'TRAV12-3*04', 'TRAV4D-4*04', 'TRAV13-4/DV7*03', 'TRAV16*01', 'TRAV9N-3*01', 'TRAV9N-2*01', 'TRAV9D-4*03', 'TRAV9D-3*01', 'TRAV9-1*01', 'TRAV9D-2*01', 'TRAV9D-4*04', 'TRAV6-1*01', 'TRAV14-1*01', 'TRAV7-2*01', 'TRAV12D-1*01', 'TRAV12-2*01', 'TRAV12N-3*01', 'TRAV16D/DV11*03', 'TRAV12-1*01', 'TRAV12-3*01', 'TRAV12D-3*02', 'TRAV12D-2*01', 'TRAV8-2*01', 'TRAV13D-2*01', 'TRAV13-4/DV7*02', 'TRAV19*01', 'TRAV21/DV12*01', 'TRAV7-5*01', 'TRAV7-5*03', 'TRAV11*01', 'TRAV11N*01', 'TRAV11*02', 'TRAV12D-1*05', 'TRAV4D-4*02', 'TRAV7-6*02', 'TRAV7-3*01', 'TRAV6D-6*03', 'TRAV7-3*02', 'TRAV15D-2/DV6D-2*04', 'TRAV15-2/DV6-2*02', 'TRAV15-2/DV6-2*01', 'TRAV15D-2/DV6D-2*01', 'TRAV15D-2/DV6D-2*02', 'TRAV15D-1/DV6D-1*04', 'TRAV15-1/DV6-1*01', 'TRAV17*02', 'TRAV17*01', 'TRAV12-1*05', 'TRAV7N-5*01', 'TRAV15D-1/DV6D-1*01', 'TRAV7-4*01', 'TRAV12D-2*02', 'TRAV9-4*01', 'TRAV9-3*01', 'TRAV6-2*01', 'TRAV8-1*01', 'TRAV13D-1*03', 'TRAV13-2*01', 'TRAV6-7/DV9*08', 'TRAV12D-1*03', 'TRAV12-3*02', 'TRAV7-6*01', 'TRAV4-2*02', 'TRAV1*01', 'TRAV5-1*01', 'TRAV8-1*02', 'TRAV15D-1/DV6D-1*05', 'TRAV7-1*01', 'TRAV9-2*01', 'TRAV13-5*01', 'TRAV6-6*03', 'TRAV12-1*02', 'TRAV12D-3*01', 'TRAV6-4*02', 'TRAV6-1*02', 'TRAV20*01', 'TRAV14D-3/DV8*05', 'TRAV18*01', 'TRAV7D-2*03', 'TRAV16*04', 'TRAV7-5*02', 'TRAV4D-2*01', 'TRAV13D-2*02', 'TRAV5D-2*01', 'TRAV5-2*01', 'TRAV16*03']) 
    assert isinstance(ng_tcrs['B'], dict) 
    assert set(ng_tcrs['B'].keys()) == set(['TRBV29*01', 'TRBV16*01', 'TRBV13-2*01', 'TRBV13-1*01', 'TRBV3*01', 'TRBV3*02', 'TRBV17*01', 'TRBV15*01', 'TRBV19*01', 'TRBV19*03', 'TRBV4*01', 'TRBV4*02', 'TRBV23*01', 'TRBV13-2*05', 'TRBV13-3*01', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV2*01', 'TRBV21*01', 'TRBV5*01', 'TRBV8*01', 'TRBV1*01', 'TRBV1*02', 'TRBV26*02', 'TRBV24*01', 'TRBV26*01', 'TRBV14*01', 'TRBV24*02', 'TRBV31*01', 'TRBV5*05', 'TRBV9*01', 'TRBV20*01', 'TRBV30*01'])  
    

def test_ng_tcrs_AOnly():
    df = clone_df_subset.iloc[0:20, :].copy()
    db = dist_b_subset.iloc[0:20, 0:20]
    da = dist_a_subset.iloc[0:20, 0:20]
    ts = TCRsubset(clone_df = df,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["A"],
            dist_a = da,
            dist_b = None)
    tm = TCRMotif(  clones_df = ts.tcr_motif_clones_df(), 
                    organism = ts.organism, 
                    chains = ts.chains, 
                    epitopes = ts.epitopes,
                    db_file = "alphabeta_db.tsv")
    ng_tcrs = dict()
    # Default Behavior
    for chain in tm.chains: #['A','B']:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref

    assert isinstance(ng_tcrs['A'], dict) #== set(['TRGV2*01', 'TRGV5P*01', 'TRGV5P*02', 'TRGV8*01', 'TRGV4*01', 'TRGV11*01', 'TRGV10*01', 'TRGV1*01', 'TRGV3*01', 'TRGV5*01', 'TRGV9*02', 'TRGV9*01'])   
    assert set(ng_tcrs['A'].keys()) == set(['TRAV14D-3/DV8*02', 'TRAV14-1*02', 'TRAV14-2*01', 'TRAV14D-1*01', 'TRAV14D-2*01', 'TRAV14D-2*02', 'TRAV14-3*01', 'TRAV14-3*02', 'TRAV14D-3/DV8*01', 'TRAV6-7/DV9*02', 'TRAV6-6*01', 'TRAV6-7/DV9*01', 'TRAV6-7/DV9*04', 'TRAV6D-6*01', 'TRAV6-5*01', 'TRAV6D-7*01', 'TRAV13-1*01', 'TRAV13D-4*01', 'TRAV13-3*01', 'TRAV10*01', 'TRAV2*01', 'TRAV6-4*01', 'TRAV3D-3*02', 'TRAV3-1*01', 'TRAV3-4*01', 'TRAV4D-4*03', 'TRAV4D-3*01', 'TRAV4D-4*01', 'TRAV4-3*01', 'TRAV4-2*01', 'TRAV4-4/DV10*01', 'TRAV4D-3*03', 'TRAV5-4*01', 'TRAV5D-4*02', 'TRAV3-3*01', 'TRAV12-3*04', 'TRAV4D-4*04', 'TRAV13-4/DV7*03', 'TRAV16*01', 'TRAV9N-3*01', 'TRAV9N-2*01', 'TRAV9D-4*03', 'TRAV9D-3*01', 'TRAV9-1*01', 'TRAV9D-2*01', 'TRAV9D-4*04', 'TRAV6-1*01', 'TRAV14-1*01', 'TRAV7-2*01', 'TRAV12D-1*01', 'TRAV12-2*01', 'TRAV12N-3*01', 'TRAV16D/DV11*03', 'TRAV12-1*01', 'TRAV12-3*01', 'TRAV12D-3*02', 'TRAV12D-2*01', 'TRAV8-2*01', 'TRAV13D-2*01', 'TRAV13-4/DV7*02', 'TRAV19*01', 'TRAV21/DV12*01', 'TRAV7-5*01', 'TRAV7-5*03', 'TRAV11*01', 'TRAV11N*01', 'TRAV11*02', 'TRAV12D-1*05', 'TRAV4D-4*02', 'TRAV7-6*02', 'TRAV7-3*01', 'TRAV6D-6*03', 'TRAV7-3*02', 'TRAV15D-2/DV6D-2*04', 'TRAV15-2/DV6-2*02', 'TRAV15-2/DV6-2*01', 'TRAV15D-2/DV6D-2*01', 'TRAV15D-2/DV6D-2*02', 'TRAV15D-1/DV6D-1*04', 'TRAV15-1/DV6-1*01', 'TRAV17*02', 'TRAV17*01', 'TRAV12-1*05', 'TRAV7N-5*01', 'TRAV15D-1/DV6D-1*01', 'TRAV7-4*01', 'TRAV12D-2*02', 'TRAV9-4*01', 'TRAV9-3*01', 'TRAV6-2*01', 'TRAV8-1*01', 'TRAV13D-1*03', 'TRAV13-2*01', 'TRAV6-7/DV9*08', 'TRAV12D-1*03', 'TRAV12-3*02', 'TRAV7-6*01', 'TRAV4-2*02', 'TRAV1*01', 'TRAV5-1*01', 'TRAV8-1*02', 'TRAV15D-1/DV6D-1*05', 'TRAV7-1*01', 'TRAV9-2*01', 'TRAV13-5*01', 'TRAV6-6*03', 'TRAV12-1*02', 'TRAV12D-3*01', 'TRAV6-4*02', 'TRAV6-1*02', 'TRAV20*01', 'TRAV14D-3/DV8*05', 'TRAV18*01', 'TRAV7D-2*03', 'TRAV16*04', 'TRAV7-5*02', 'TRAV4D-2*01', 'TRAV13D-2*02', 'TRAV5D-2*01', 'TRAV5-2*01', 'TRAV16*03']) 

def test_ng_tcrs_BOnly():
    df = clone_df_subset.iloc[0:20, :].copy()
    db = dist_b_subset.iloc[0:20, 0:20]
    da = dist_a_subset.iloc[0:20, 0:20]
    ts = TCRsubset(clone_df = df,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["B"],
            dist_a = da,
            dist_b = db)
    tm = TCRMotif(  clones_df = ts.tcr_motif_clones_df(), 
                    organism = ts.organism, 
                    chains = ts.chains, 
                    epitopes = ts.epitopes,
                    db_file = "alphabeta_db.tsv")
    ng_tcrs = dict()
    # Default Behavior
    for chain in tm.chains: #['A','B']:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref

    assert isinstance(ng_tcrs['B'], dict) 
    assert set(ng_tcrs['B'].keys()) == set(['TRBV29*01', 'TRBV16*01', 'TRBV13-2*01', 'TRBV13-1*01', 'TRBV3*01', 'TRBV3*02', 'TRBV17*01', 'TRBV15*01', 'TRBV19*01', 'TRBV19*03', 'TRBV4*01', 'TRBV4*02', 'TRBV23*01', 'TRBV13-2*05', 'TRBV13-3*01', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV2*01', 'TRBV21*01', 'TRBV5*01', 'TRBV8*01', 'TRBV1*01', 'TRBV1*02', 'TRBV26*02', 'TRBV24*01', 'TRBV26*01', 'TRBV14*01', 'TRBV24*02', 'TRBV31*01', 'TRBV5*05', 'TRBV9*01', 'TRBV20*01', 'TRBV30*01'])  
    

def test_ng_tcrs_GD():
    ts = TCRsubset(clone_df = clone_df_subset_d,            
            organism = "human",
            epitopes = ["X"] ,
            epitope = "X",
            chains = ["delta"],
            dist_d = dist_d_subset)

    tcr_motif_delta_input = pd.DataFrame({'subject': {53: 'SRR5130260.1',  54: 'SRR5130260.1',  55: 'SRR5130260.1',  56: 'SRR5130260.1',  57: 'SRR5130260.1',  58: 'SRR5130260.1',  59: 'SRR5130260.1',  60: 'SRR5130260.1',  61: 'SRR5130260.1',  62: 'SRR5130260.1',  63: 'SRR5130260.1',  64: 'SRR5130260.1',  65: 'SRR5130260.1'}, 'epitope': {53: 'X',  54: 'X',  55: 'X',  56: 'X',  57: 'X',  58: 'X',  59: 'X',  60: 'X',  61: 'X',  62: 'X',  63: 'X',  64: 'X',  65: 'X'}, 'va_rep': {53: 'TRGV1*01',  54: 'TRGV1*01',  55: 'TRGV1*01',  56: 'TRGV1*01',  57: 'TRGV1*01',  58: 'TRGV1*01',  59: 'TRGV1*01',  60: 'TRGV1*01',  61: 'TRGV1*01',  62: 'TRGV1*01',  63: 'TRGV1*01',  64: 'TRGV1*01',  65: 'TRGV1*01'}, 'ja_rep': {53: 'TRGJ1*01',  54: 'TRGJ1*01',  55: 'TRGJ1*01',  56: 'TRGJ1*01',  57: 'TRGJ1*01',  58: 'TRGJ1*01',  59: 'TRGJ1*01',  60: 'TRGJ1*01',  61: 'TRGJ1*01',  62: 'TRGJ1*01',  63: 'TRGJ1*01',  64: 'TRGJ1*01',  65: 'TRGJ1*01'}, 'vb_rep': {53: 'TRDV2*01',  54: 'TRDV2*01',  55: 'TRDV2*01',  56: 'TRDV2*01',  57: 'TRDV2*01',  58: 'TRDV2*01',  59: 'TRDV2*01',  60: 'TRDV2*01',  61: 'TRDV2*01',  62: 'TRDV2*01',  63: 'TRDV2*01',  64: 'TRDV2*01',  65: 'TRDV2*01'}, 'jb_rep': {53: 'TRDJ1*01',  54: 'TRDJ1*01',  55: 'TRDJ1*01',  56: 'TRDJ1*01',  57: 'TRDJ2*01',  58: 'TRDJ2*01',  59: 'TRDJ2*01',  60: 'TRDJ2*01',  61: 'TRDJ2*01',  62: 'TRDJ2*01',  63: 'TRDJ2*01',  64: 'TRDJ2*01',  65: 'TRDJ2*01'}, 'cdr3a': {53: 'CATWAKNYYKKLF',  54: 'CATWAKNYYKKLF',  55: 'CATWAKNYYKKLF',  56: 'CATWAKNYYKKLF',  57: 'CATWAKNYYKKLF',  58: 'CATWAKNYYKKLF',  59: 'CATWAKNYYKKLF',  60: 'CATWAKNYYKKLF',  61: 'CATWAKNYYKKLF',  62: 'CATWAKNYYKKLF',  63: 'CATWAKNYYKKLF',  64: 'CATWAKNYYKKLF',  65: 'CATWAKNYYKKLF'}, 'cdr3b': {53: 'CACHRGTDTDKLIF',  54: 'CACDKNGGYVRYTDKLIF',  55: 'CACDTVGIPDKLIF',  56: 'CACVRLPLRGRPYTDKLIF',  57: 'CACDNWGALTAQLFF',  58: 'CACDTILGDITLTAQLFF',  59: 'CACDTGRGTLTAQLFF',  60: 'CACDTWGMTAQLFF',  61: 'CACDTGGALTAQLFF',  62: 'CACDIRDTRVLTAQLFF',  63: 'CACDIVLGDPSLTAQLFF',  64: 'CACDHLLGDTAQLFF',  65: 'CACDPVTGGSLTAQLFF'}})
    assert np.all( ts.tcr_motif_clones_df(gdmode=True)   == tcr_motif_delta_input)

    tm = TCRMotif(  clones_df = ts.tcr_motif_clones_df(gdmode = True), 
                    organism = ts.organism, 
                    chains = ts.chains, 
                    epitopes = ts.epitopes,
                    db_file = "gammadelta_db.tsv")
    
    ng_tcrs = dict()
    for chain in tm.chains:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref
    
    # Notice that 'B' represents 'delta'
    # Notice that 'A' represents 'gamma'

    ng_tcrs = dict()
    for chain in ['A','B']:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref

    assert set(ng_tcrs['B'].keys()) == set(['TRDV2*01', 'TRDV3*01', 'TRDV1*01', 'TRAV38-2/DV8*01', 'TRAV22*01', 'TRAV29/DV5*01', 'TRAV41*01', 'TRAV39*01', 'TRAV14/DV4*02', 'TRAV40*01', 'TRAV23/DV6*01', 'TRAV34*01', 'TRAV38-1*01', 'TRAV26-2*01', 'TRAV26-1*01', 'TRAV19*01', 'TRAV35*01', 'TRAV17*01', 'TRAV20*01', 'TRAV21*01', 'TRAV36/DV7*01', 'TRAV14/DV4*01', 'TRAV9-2*01', 'TRAV24*01', 'TRAV30*01', 'TRAV38-1*03', 'TRAV38-1*02', 'TRAV36/DV7*02', 'TRAV12-3*01', 'TRAV27*01', 'TRAV8-3*01', 'TRAV16*01', 'TRAV13-1*01', 'TRAV30*03', 'TRAV8-4*01', 'TRAV8-4*06', 'TRAV8-2*01', 'TRAV10*01', 'TRAV8-4*07', 'TRAV8-6*01', 'TRAV25*01', 'TRAV9-1*01', 'TRAV36/DV7*03', 'TRAV6*01', 'TRAV13-2*01', 'TRAV8-7*01'])  
    assert set(ng_tcrs['A'].keys()) == set(['TRGV2*01', 'TRGV5P*01', 'TRGV5P*02', 'TRGV8*01', 'TRGV4*01', 'TRGV11*01', 'TRGV10*01', 'TRGV1*01', 'TRGV3*01', 'TRGV5*01', 'TRGV9*02', 'TRGV9*01'])   



def test_ng_tcrs_DOnly():
    ts = TCRsubset(clone_df = clone_df_subset_d,            
            organism = "human",
            epitopes = ["X"] ,
            epitope = "X",
            chains = ["delta"],
            dist_d = dist_d_subset)

    tcr_motif_delta_input = pd.DataFrame({'subject': {53: 'SRR5130260.1',  54: 'SRR5130260.1',  55: 'SRR5130260.1',  56: 'SRR5130260.1',  57: 'SRR5130260.1',  58: 'SRR5130260.1',  59: 'SRR5130260.1',  60: 'SRR5130260.1',  61: 'SRR5130260.1',  62: 'SRR5130260.1',  63: 'SRR5130260.1',  64: 'SRR5130260.1',  65: 'SRR5130260.1'}, 'epitope': {53: 'X',  54: 'X',  55: 'X',  56: 'X',  57: 'X',  58: 'X',  59: 'X',  60: 'X',  61: 'X',  62: 'X',  63: 'X',  64: 'X',  65: 'X'}, 'va_rep': {53: 'TRGV1*01',  54: 'TRGV1*01',  55: 'TRGV1*01',  56: 'TRGV1*01',  57: 'TRGV1*01',  58: 'TRGV1*01',  59: 'TRGV1*01',  60: 'TRGV1*01',  61: 'TRGV1*01',  62: 'TRGV1*01',  63: 'TRGV1*01',  64: 'TRGV1*01',  65: 'TRGV1*01'}, 'ja_rep': {53: 'TRGJ1*01',  54: 'TRGJ1*01',  55: 'TRGJ1*01',  56: 'TRGJ1*01',  57: 'TRGJ1*01',  58: 'TRGJ1*01',  59: 'TRGJ1*01',  60: 'TRGJ1*01',  61: 'TRGJ1*01',  62: 'TRGJ1*01',  63: 'TRGJ1*01',  64: 'TRGJ1*01',  65: 'TRGJ1*01'}, 'vb_rep': {53: 'TRDV2*01',  54: 'TRDV2*01',  55: 'TRDV2*01',  56: 'TRDV2*01',  57: 'TRDV2*01',  58: 'TRDV2*01',  59: 'TRDV2*01',  60: 'TRDV2*01',  61: 'TRDV2*01',  62: 'TRDV2*01',  63: 'TRDV2*01',  64: 'TRDV2*01',  65: 'TRDV2*01'}, 'jb_rep': {53: 'TRDJ1*01',  54: 'TRDJ1*01',  55: 'TRDJ1*01',  56: 'TRDJ1*01',  57: 'TRDJ2*01',  58: 'TRDJ2*01',  59: 'TRDJ2*01',  60: 'TRDJ2*01',  61: 'TRDJ2*01',  62: 'TRDJ2*01',  63: 'TRDJ2*01',  64: 'TRDJ2*01',  65: 'TRDJ2*01'}, 'cdr3a': {53: 'CATWAKNYYKKLF',  54: 'CATWAKNYYKKLF',  55: 'CATWAKNYYKKLF',  56: 'CATWAKNYYKKLF',  57: 'CATWAKNYYKKLF',  58: 'CATWAKNYYKKLF',  59: 'CATWAKNYYKKLF',  60: 'CATWAKNYYKKLF',  61: 'CATWAKNYYKKLF',  62: 'CATWAKNYYKKLF',  63: 'CATWAKNYYKKLF',  64: 'CATWAKNYYKKLF',  65: 'CATWAKNYYKKLF'}, 'cdr3b': {53: 'CACHRGTDTDKLIF',  54: 'CACDKNGGYVRYTDKLIF',  55: 'CACDTVGIPDKLIF',  56: 'CACVRLPLRGRPYTDKLIF',  57: 'CACDNWGALTAQLFF',  58: 'CACDTILGDITLTAQLFF',  59: 'CACDTGRGTLTAQLFF',  60: 'CACDTWGMTAQLFF',  61: 'CACDTGGALTAQLFF',  62: 'CACDIRDTRVLTAQLFF',  63: 'CACDIVLGDPSLTAQLFF',  64: 'CACDHLLGDTAQLFF',  65: 'CACDPVTGGSLTAQLFF'}})
    assert np.all( ts.tcr_motif_clones_df(gdmode=True)   == tcr_motif_delta_input)

    tm = TCRMotif(  clones_df = ts.tcr_motif_clones_df(gdmode = True), 
                    organism = ts.organism, 
                    chains = ts.chains, 
                    epitopes = ts.epitopes,
                    db_file = "gammadelta_db.tsv")
    
    ng_tcrs = dict()
    for chain in tm.chains:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref
    
    # Notice that 'B' represents 'delta'
    # Notice that 'A' represents 'gamma'

    ng_tcrs = dict()
    for chain in tm.chains:
        next_gen_ref = tm.generate_background_set(chain = chain,
                            ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
                            ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
        ng_tcrs[chain] = next_gen_ref

    assert set(ng_tcrs['B'].keys()) == set(['TRDV2*01', 'TRDV3*01', 'TRDV1*01', 'TRAV38-2/DV8*01', 'TRAV22*01', 'TRAV29/DV5*01', 'TRAV41*01', 'TRAV39*01', 'TRAV14/DV4*02', 'TRAV40*01', 'TRAV23/DV6*01', 'TRAV34*01', 'TRAV38-1*01', 'TRAV26-2*01', 'TRAV26-1*01', 'TRAV19*01', 'TRAV35*01', 'TRAV17*01', 'TRAV20*01', 'TRAV21*01', 'TRAV36/DV7*01', 'TRAV14/DV4*01', 'TRAV9-2*01', 'TRAV24*01', 'TRAV30*01', 'TRAV38-1*03', 'TRAV38-1*02', 'TRAV36/DV7*02', 'TRAV12-3*01', 'TRAV27*01', 'TRAV8-3*01', 'TRAV16*01', 'TRAV13-1*01', 'TRAV30*03', 'TRAV8-4*01', 'TRAV8-4*06', 'TRAV8-2*01', 'TRAV10*01', 'TRAV8-4*07', 'TRAV8-6*01', 'TRAV25*01', 'TRAV9-1*01', 'TRAV36/DV7*03', 'TRAV6*01', 'TRAV13-2*01', 'TRAV8-7*01'])  


def test_non_default_manual_generation_of_ng_tcrsAB():
    """ Notices the usage of non default mode and specification of files"""
    df = clone_df_subset.iloc[0:20, :].copy()
    db = dist_b_subset.iloc[0:20, 0:20]
    da = dist_a_subset.iloc[0:20, 0:20]
    ts=TCRsubset(clone_df = df,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["A","B"],
            dist_a = da,
            dist_b = db,
            default_mode = False)
    ts.ng_tcrs['B'] = ts.generate_background_set(chain = ['B'],ng_log_path = 'tcrdist/db/alphabeta_db.tsv_files',ng_log_file = 'new_nextgen_chains_mouse_B.tsv')
    ts.ng_tcrs['A'] = ts.generate_background_set(chain = ['A'],ng_log_path = 'tcrdist/db/alphabeta_db.tsv_files',ng_log_file = 'new_nextgen_chains_mouse_A.tsv')    
    
    motif_df = ts.find_motif()
    assert isinstance(motif_df, pd.DataFrame)
    assert isinstance(ts.motif_df, pd.DataFrame)
    assert motif_df.shape[0] > 1









# @pytest.mark.skip()
# def develoment_test_only():
#     import pytest
#     import pandas as pd
#     from tcrdist.subset import TCRsubset
#     from tcrdist.tests.my_test_subset import dist_d_subset, clone_df_subset_d
#     from tcrdist.tests.my_test_subset import dist_a_subset, dist_b_subset, clone_df_subset 
#     import numpy as np
#     from tcrdist.cdr3_motif import TCRMotif
#     from tcrdist import paths 

#     assert isinstance(dist_d_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset_d, pd.DataFrame)
    
#     clone_df_subset_d['epitope'] = "X"
    
#     ts = TCRsubset(clone_df = clone_df_subset_d,            
#             organism = "human",
#             epitopes = ["X"] ,
#             epitope = "X",
#             chains = ["delta"],
#             dist_d = dist_d_subset)

#     [(x, x in ts.tcr_motif_clones_df(gdmode = True).columns) for x in ['epitope','va_rep','ja_rep','vb_rep','jb_rep','cdr3a','cdr3b']]

#     tm = TCRMotif(clones_df = ts.tcr_motif_clones_df(gdmode = True), 
#                     organism = ts.organism, 
#                     chains = ts.chains, 
#                     epitopes = ts.epitopes,
#                     db_file = "gammadelta_db.tsv")  
    
#     # default_behaivor
#     ng_tcrs = dict()
#     for chain in tm.chains:
#         next_gen_ref = tm.generate_background_set(chain = chain,
#                             ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
#                             ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,chain) )
#         # next_gen_a = tm.generate_background_set(chain = 'A',
#         #                     ng_log_path = paths.path_to_current_db_files(db_file = tm.db_file),
#         #                     ng_log_file = 'new_nextgen_chains_{}_{}.tsv'.format(tm.organism,'A') )
#         ng_tcrs[chain] = next_gen_ref

#     #ng_tcrs = {**next_gen_b, **next_gen_a}
    
#     #assert set(ng_tcrs['A'].keys()) == set(['TRGV2*01', 'TRGV5P*01', 'TRGV5P*02', 'TRGV8*01', 'TRGV4*01', 'TRGV11*01', 'TRGV10*01', 'TRGV1*01', 'TRGV3*01', 'TRGV5*01', 'TRGV9*02', 'TRGV9*01'])    
#     assert set(ng_tcrs['B'].keys()) == set(['TRDV2*01', 'TRDV3*01', 'TRDV1*01', 'TRAV38-2/DV8*01', 'TRAV22*01', 'TRAV29/DV5*01', 'TRAV41*01', 'TRAV39*01', 'TRAV14/DV4*02', 'TRAV40*01', 'TRAV23/DV6*01', 'TRAV34*01', 'TRAV38-1*01', 'TRAV26-2*01', 'TRAV26-1*01', 'TRAV19*01', 'TRAV35*01', 'TRAV17*01', 'TRAV20*01', 'TRAV21*01', 'TRAV36/DV7*01', 'TRAV14/DV4*01', 'TRAV9-2*01', 'TRAV24*01', 'TRAV30*01', 'TRAV38-1*03', 'TRAV38-1*02', 'TRAV36/DV7*02', 'TRAV12-3*01', 'TRAV27*01', 'TRAV8-3*01', 'TRAV16*01', 'TRAV13-1*01', 'TRAV30*03', 'TRAV8-4*01', 'TRAV8-4*06', 'TRAV8-2*01', 'TRAV10*01', 'TRAV8-4*07', 'TRAV8-6*01', 'TRAV25*01', 'TRAV9-1*01', 'TRAV36/DV7*03', 'TRAV6*01', 'TRAV13-2*01', 'TRAV8-7*01'])  



# @pytest.mark.skip()
# def test_that_creation_of_TCRsubset_delta_only():
#     """
#     Tests creation of TCRsubset from HUMAN delta chain bulk sequences 
#     """
#     import pytest
#     import pandas as pd
#     from tcrdist.subset import TCRsubset
#     from tcrdist.tests.my_test_subset import dist_d_subset, clone_df_subset_d
#     from tcrdist.tests.my_test_subset import dist_a_subset, dist_b_subset, clone_df_subset 


#     assert isinstance(dist_d_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset_d, pd.DataFrame)
    
#     clone_df_subset_d['epitope'] = "X"
    
#     ts = TCRsubset(clone_df = clone_df_subset_d,            
#             organism = "human",
#             epitopes = ["X"] ,
#             epitope = "X",
#             chains = ["delta"],
#             dist_d = dist_d_subset)

# @pytest.mark.skip()
# def test_that_creation_of_TCRsubset_plus_find_motif_delta_only():
#     """
#     Tests createion of TCRsubset from HUMAN delta chain bulk sequences 
#     Then try
#     """
#     assert isinstance(dist_d_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset_d, pd.DataFrame)
    
#     clone_df_subset_d['epitope'] = "X"
    
#     ts = TCRsubset(clone_df = clone_df_subset_d,            
#             organism = "human",
#             epitopes = ["X"] ,
#             epitope = "X",
#             chains = ["delta"],
#             dist_d = dist_d_subset)

#     ts.find_motif()


# @pytest.mark.skip()
# def test_initialization_of_TCRsubset_alpha_beta_case_plus_motif_finding():
#     """
#     Test that we can create a TCRsubset from paired input files
#     """
#     assert isinstance(dist_a_subset, pd.DataFrame)
#     assert isinstance(dist_b_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset, pd.DataFrame)
#     df = clone_df_subset.iloc[0:20, :].copy()
#     db = dist_b_subset.iloc[0:20, 0:20]
#     da = dist_a_subset.iloc[0:20, 0:20]
#     ts=TCRsubset(clone_df = df,            
#             organism = "mouse",
#             epitopes = ["PA"] ,
#             epitope = "PA",
#             chains = ["A","B"],
#             dist_a = da,
#             dist_b = db)
#     motif_df = ts.find_motif()
#     assert isinstance(motif_df, pd.DataFrame)
#     assert isinstance(ts.motif_df, pd.DataFrame)

# @pytest.mark.skip()
# def test_initialization_of_TCRsubset_beta_case():
#     """
#     Test that we can create a TCRsubset from just beta chain files.
#     For this to work, fake alpha sequences will have to be created:
#     """
#     assert isinstance(dist_a_subset, pd.DataFrame)
#     assert isinstance(dist_b_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset, pd.DataFrame)
#     df = clone_df_subset[['clone_id', 'subject', 'epitope', 
#                         'v_b_gene', 'j_b_gene', 
#                         'cdr3_b_aa', 'cdr1_b_aa', 'cdr2_b_aa', 
#                         'pmhc_b_aa', 'cdr3_b_nucseq', 'count',
#                         'vb_countreps', 'jb_countreps','vb_gene', 'jb_gene']]
    
#     ts = TCRsubset(  clone_df = df,            
#                 organism = "mouse",
#                 epitopes = ["PA"] ,
#                 epitope = "PA",
#                 chains = ["beta"],
#                 dist_b = dist_b_subset)
    
#     assert(isinstance(ts.clone_df.va_gene, pd.Series))
#     assert(ts.clone_df.va_gene.iloc[0] == "TRAV10*01")


# @pytest.mark.skip()
# def test_that_test_dataframes_exist():
#     assert isinstance(dist_d_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset_d, pd.DataFrame)



# @pytest.mark.skip()
# def test_initialization_of_TCRsubset_alpha_beta_case_with_gamma_delta_chain_names():
#     """
#     Test that we can create a TCRsubset from paired input files
#     """
#     assert isinstance(dist_a_subset, pd.DataFrame)
#     assert isinstance(dist_b_subset, pd.DataFrame)
#     assert isinstance(clone_df_subset, pd.DataFrame)
#     TCRsubset(clone_df = clone_df_subset,            
#             organism = "mouse",
#             epitopes = ["PA"] ,
#             epitope = "PA",
#             chains = ["gamma","delta"],
#             dist_a = dist_a_subset,
#             dist_b = dist_b_subset)
        

    
