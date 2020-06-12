"""
2020-05-29 test_repertoire_ui.py includes integration tests for user facing commands from TCRrep
"""
import pytest
import numpy as np
import tcrdist as td
from tcrdist.repertoire import TCRrep
import pandas as pd
import numpy as np
import os
from os.path import join as opj
import sys
import multiprocessing
from zipdist.zip2 import Zipdist2

def run_TCRrep_func_tcrdist2(chains = ['beta','alpha'], metric = "nw"):
    cpu = multiprocessing.cpu_count()
    # really basic example
    df = pd.read_csv(opj('tcrdist', 'datasets', 'dash.csv'))
    df = df[df.epitope.isin(['NP'])]
    tr = TCRrep(cell_df=df, chains=chains, organism='mouse')
    tr.tcrdist2(processes = cpu, metric = metric, dump = True, reduce = True, save=False)
    return tr

def test_TCRrep_func_tcrdist2_save_manual_rebuild(chains = ['beta','alpha'], metric = "nw"):
    cpu = multiprocessing.cpu_count()
    # really basic example
    df = pd.read_csv(opj('tcrdist', 'datasets', 'dash.csv'))
    df = df[df.epitope.isin(['NP'])]
    tr = TCRrep(cell_df=df, chains=chains, organism='mouse')
    tr.tcrdist2(processes = cpu, metric = metric, dump = True, reduce = True, save=True, dest = "default_archive", dest_tar_name = "default_archive.tar.gz" )
    # Cleanup folder that you just made
    os.system("rm -rf myTCRrep_archive")
    # Rebuild
    tr = TCRrep(cell_df=df.iloc[0:0,:], chains=chains, organism='mouse')
    z = Zipdist2(name = "test_only", target = tr)
    z._build(dest_tar = "default_archive.tar.gz", target = tr)
    assert isinstance(tr.paired_tcrdist, np.ndarray )
    assert isinstance(tr.pw_tcrdist, np.ndarray )
    assert np.array_equal(tr.pw_tcrdist, tr.paired_tcrdist)


def test_TCRrep_func_tcrdist2_save_auto_rebuild(chains = ['beta','alpha'], metric = "nw"):
    cpu = multiprocessing.cpu_count()
    # really basic example
    df = pd.read_csv(opj('tcrdist', 'datasets', 'dash.csv'))
    df = df[df.epitope.isin(['NP'])]
    tr = TCRrep(cell_df=df, chains=chains, organism='mouse')
    tr.tcrdist2(processes = cpu, metric = metric, dump = True, reduce = True, save= True, dest = "default_archive", dest_tar_name = "default_archive.tar.gz" )
    # Cleanup folder that you just made
    os.system("rm -rf myTCRrep_archive")
    # Rebuild
    tr = TCRrep(cell_df=df.iloc[0:0,:], chains=chains, organism='mouse')
    tr.rebuild(dest_tar_name = "default_archive.tar.gz")

    tr_compare = run_TCRrep_func_tcrdist2(chains = ['beta','alpha'], metric = "nw")

    attendance = {k : k in tr.__dict__.keys() for k in tr_compare.__dict__.keys()}
    assert attendance['pw_tcrdist']
    assert attendance['pw_alpha']
    assert attendance['pw_beta']
    assert attendance['clone_df']
    assert attendance['cell_df']
    # only compare things in common
    shared_attributes = [k for k in attendance.keys() if attendance[k]]
    for k in shared_attributes:
        if not isinstance(getattr(tr, k), pd.DataFrame):
            if not isinstance(getattr(tr, k), dict):
                assert np.all(getattr(tr, k) == getattr(tr_compare, k))
            else:
                assert set(getattr(tr, k).keys()) ==  set(getattr(tr_compare, k).keys())

    assert set(tr.all_genes['mouse'].keys()) == set(tr_compare.all_genes['mouse'].keys())
    assert set(tr.all_genes['human'].keys()) == set(tr_compare.all_genes['human'].keys())


def test_TCRrep_func_trcrdist2_AB_nw():
    tr = run_TCRrep_func_tcrdist2(['beta','alpha'])
    assert tr.pw_beta.shape[0] == tr.pw_beta.shape[1]
    assert tr.pw_alpha.shape[0] == tr.pw_alpha.shape[1]
    assert tr.pw_beta.shape[0] == tr.pw_alpha.shape[1]
    assert isinstance(tr.pw_alpha, np.ndarray)
    assert isinstance(tr.pw_beta, np.ndarray)
    assert isinstance(tr.paired_tcrdist, np.ndarray)
    assert np.array_equal(tr.pw_tcrdist, tr.paired_tcrdist)
    #assert tr.paired_tcrdist is tr.pw_tcrdist


def test_TCRrep_func_trcrdist2_AB_nwhamming():
    tr = run_TCRrep_func_tcrdist2(['beta','alpha'], metric = "hamming")
    assert tr.pw_beta.shape[0] == tr.pw_beta.shape[1]
    assert tr.pw_alpha.shape[0] == tr.pw_alpha.shape[1]
    assert tr.pw_beta.shape[0] == tr.pw_alpha.shape[1]
    assert isinstance(tr.pw_alpha, np.ndarray)
    assert isinstance(tr.pw_beta, np.ndarray)
    assert isinstance(tr.paired_tcrdist, np.ndarray)
    assert isinstance(tr.pw_tcrdist, np.ndarray)
    #assert tr.paired_tcrdist is tr.pw_tcrdist


def test_TCRrep_func_trcrdist2_B_nw():
    tr = run_TCRrep_func_tcrdist2(['beta'])
    assert tr.pw_beta.shape[0] == tr.pw_beta.shape[1]
    with pytest.raises(AttributeError):
        tr.pw_alpha is None
    assert isinstance(tr.pw_beta, np.ndarray)
    assert isinstance(tr.paired_tcrdist, np.ndarray)
    assert isinstance(tr.pw_tcrdist, np.ndarray)
    #assert tr.paired_tcrdist is tr.pw_tcrdist #NO LONGER TRUE AFTER reduce


def test_TCRrep_func_trcrdist2_A_nw():
    tr = run_TCRrep_func_tcrdist2(['alpha'])
    assert tr.pw_alpha.shape[0] == tr.pw_alpha.shape[1]
    assert isinstance(tr.pw_alpha, np.ndarray)
    with pytest.raises(AttributeError):
        tr.pw_beta is None    
    assert isinstance(tr.paired_tcrdist, np.ndarray)
    assert isinstance(tr.pw_tcrdist, np.ndarray)
    #assert tr.paired_tcrdist is tr.pw_tcrdist
    #assert np.all(tr.paired_tcrdist == tr.pw_alpha)
    
    