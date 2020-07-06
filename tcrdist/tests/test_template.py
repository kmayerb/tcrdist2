import os
import pandas as pd
import pytest
import numpy as np

from tcrdist import cluster
from tcrsampler.sampler import TCRsampler
import collections

ts = TCRsampler(default_background = "wirasinha_mouse_beta_t_1.tsv.sampler.tsv")

fn = os.path.join('tcrdist','tests','dash_pa_clone_df.csv')
test_clone_df = pd.read_csv(fn)

fn = os.path.join('tcrdist','tests','dash_pa_cluster_df.csv')
test_cluster_df = pd.read_csv(fn)

def test_that_inputs_exist():
    assert isinstance(test_clone_df, pd.DataFrame)
    assert isinstance(test_cluster_df, pd.DataFrame)

def test_df_to_toc_html():
    cluster.dataframe_to_toc_html(df = test_cluster_df)

def test_cluster_df_to_cluster_html():
    motifs_list = list()
    for i,row in test_cluster_df.head(10).iterrows():
        cl_attrs = cluster._get_cluster_attributes(row = row, df= test_clone_df, )
        motifs = cluster._get_palmotifs(cluster_attributes = cl_attrs, 
                        sampler = ts, 
                        write = True, 
                        dest = 'static')
        motifs_list.append(motifs)
        cluster.cluster_dataframe_to_cluster_html(df = cl_attrs.ns_df, 
                                                cluster_id = cl_attrs.cluster_id,
                                                cols = ["cdr3_b_aa", 'v_b_gene', 'j_b_gene', 'subject','epitope'])
    
    
    # Write the table of contents
    df = pd.DataFrame({'cluster_id' : [x.cluster_attributes.cluster_id for x in motifs_list],
                       'centroid' : [x.cluster_attributes.centroid for x in motifs_list],
                       'loglik': [-1*np.sum(x.pal_stat.loglik) for x in motifs_list],
                       'warning': [x.warning for x in motifs_list] })
    
    dfx = test_cluster_df.merge(df, how = 'left', on = 'cluster_id')

    cluster.dataframe_to_toc_html(df = dfx, cols=['cluster_id','neighbors','K_neighbors','centroid','loglik','warning'] )
