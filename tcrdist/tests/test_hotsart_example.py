import pytest
   # All code in one block minus, except plotting for the purposes of testing.
def test_hot_start_example_in_full():
    """
    This is the code that makes up the HotStart example in the docs
    """
    # basic imports
    import os
    import pandas as pd
    import numpy as np
    #import IPython

    # tcrdist classes
    from tcrdist.repertoire import TCRrep
    from tcrdist.subset import TCRsubset
    from tcrdist.cdr3_motif import TCRMotif
    from tcrdist.storage import StoreIOMotif, StoreIOEntropy

    # tcrdist functions
    from tcrdist import plotting
    from tcrdist.mappers import populate_legacy_fields

    # scipy functions for clustering
    from scipy.spatial import distance    
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster 

    # sklearn functions for low-dimensional embeddings 
    from sklearn.manifold import TSNE, MDS

    # plotnine to allow grammar of graphics plotting akin to R's ggplot2
    #import plotnine as gg

    #1 load data, subset to receptors recognizing "PA" epitope
    tcrdist2_df = pd.read_csv(os.path.join("tcrdist","test_files_compact","dash.csv"))
    tcrdist2_df = tcrdist2_df[tcrdist2_df.epitope == "PA"].copy()

    #2 create instance of TCRrep class, initializes input as tr.cell_df attribute
    tr = TCRrep(cell_df = tcrdist2_df, chains = ['alpha','beta'],organism = "mouse")

    #3 Infer CDR1,CDR2,CDR2.5 (a.k.a. phmc) from germline v-genes
    tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
    tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)

    #4 Define index columns for determining unique clones.
    tr.index_cols = ['clone_id', 'subject', 'epitope',
                    'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                    'cdr3_a_aa', 'cdr3_b_aa',
                    'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                    'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                    'cdr3_b_nucseq', 'cdr3_a_nucseq']

    #4 Deduplicate based on index cols, creating tr.clone_df attribute
    tr.deduplicate()

    #5 calculate tcrdists by method in Dash et al. 
    tr._tcrdist_legacy_method_alpha_beta()

    #6 Check that sum of alpah-chain and beta-chain distance matrices equal paired_tcrdist
    distA = tr.dist_a
    distB = tr.dist_b
    assert np.all(((distA + distB) - tr.paired_tcrdist) == 0)


    # Cluster
    from scipy.spatial import distance    
    from scipy.cluster.hierarchy import linkage, dendrogram, fcluster 
    compressed_dmat = distance.squareform(tr.paired_tcrdist, force = "vector")
    Z = linkage(compressed_dmat, method = "complete")
    den = dendrogram(Z, color_threshold = np.inf, no_plot = True)
    cluster_index = fcluster(Z, t = 20, criterion = "maxclust")
    assert len(cluster_index) == tr.clone_df.shape[0]
    assert len(cluster_index) == tr.paired_tcrdist.shape[0]
    tr.clone_df['cluster_index'] = cluster_index

    # Subset to Cluster 5
    criteria = (cluster_index == 5) 
    clone_df_subset = tr.clone_df[criteria]
    clone_df_subset = clone_df_subset[clone_df_subset.epitope == "PA"].copy()
    dist_a_subset = tr.dist_a.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()
    dist_b_subset = tr.dist_b.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()

    clone_df_subset = populate_legacy_fields(df = clone_df_subset, chains =['alpha', 'beta'])

    ts = TCRsubset(clone_df_subset,
                organism = "mouse",
                epitopes = ["PA"] ,
                epitope = "PA",
                chains = ["A","B"],
                dist_a = dist_a_subset,
                dist_b = dist_b_subset)

    # Find Motifs 
    if os.path.isfile(os.path.join("tcrdist", "test_files_compact","dash_PA_cluster_5_motifs.csv")):
        ts.motif_df = pd.read_csv(os.path.join("tcrdist", "test_files_compact","dash_PA_cluster_5_motifs.csv"))
    else:
        motif_df = ts.find_motif()

    # Save Motifs
    ts.motif_df.to_csv(os.path.join("tcrdist", "test_files_compact","dash_PA_cluster_5_motifs.csv"), index = False)

    # Preprocess Motifs 
    motif_list_a = list()
    motif_logos_a = list()
    for i,row in ts.motif_df[ts.motif_df.ab == "A"].iterrows():
        StoreIOMotif_instance = ts.eval_motif(row)
        motif_list_a.append(StoreIOMotif_instance)
        motif_logos_a.append(plotting.plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600))
            
    motif_list_b = list()
    motif_logos_b = list()
    for i,row in ts.motif_df[ts.motif_df.ab == "B"].iterrows():
        StoreIOMotif_instance = ts.eval_motif(row)
        motif_list_b.append(StoreIOMotif_instance)
        motif_logos_b.append(plotting.plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600))


