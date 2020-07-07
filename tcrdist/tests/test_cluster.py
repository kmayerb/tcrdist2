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
