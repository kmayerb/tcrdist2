import pytest

def test_CompleteExample_with_TCRMotif_Invoked_From_within_TCRsubset():
    import pandas as pd
    import numpy as np
    import tcrdist as td
    #import IPython

    from tcrdist import mappers
    from tcrdist.repertoire import TCRrep
    from tcrdist.cdr3_motif import TCRMotif
    from tcrdist.subset import TCRsubset
    from tcrdist.storage import StoreIOMotif, StoreIOEntropy
    from tcrdist.plotting import plot_pwm

    tcrdist_clone_fn = 'tcrdist/test_files/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'
    tcrdist_clone_df = pd.read_csv(tcrdist_clone_fn, sep = "\t")               #1

    ind = (tcrdist_clone_df.epitope == "PA") | (tcrdist_clone_df.epitope == "F2")
    tcrdist_clone_df = tcrdist_clone_df[ind].copy()

    mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping                     #3
    tcrdist2_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df,         #4
                                                mapping = mapping)



    #1
    tr = TCRrep(cell_df = tcrdist2_df, organism = "mouse")

    #2
    tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
    tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)

    #3
    tr.index_cols = ['clone_id', 'subject', 'epitope',
                     'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                     'cdr3_a_aa', 'cdr3_b_aa',
                     'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                     'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                     'cdr3_b_nucseq', 'cdr3_a_nucseq',
                     'va_countreps', 'ja_countreps',
                     'vb_countreps', 'jb_countreps',
                     'va_gene', 'vb_gene',
                     'ja_gene', 'jb_gene']

    #4
    tr.deduplicate()

    #5
    tr._tcrdist_legacy_method_alpha_beta()

    #6
    distA = tr.dist_a
    distB = tr.dist_b
    assert np.all(((distA + distB) - tr.paired_tcrdist) == 0)

    # 1
    criteria = tr.clone_df.epitope == "PA"
    clone_df_subset = tr.clone_df[criteria]

    # 2
    distA_subset = distA.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()
    distB_subset = distB.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()

    # 3
    ts = TCRsubset(clone_df_subset,
                 organism = "mouse",
                 epitopes = ["PA"],
                 epitope = "PA",
                 chains = ["A","B"],
                 dist_a = distA_subset,
                 dist_b = distB_subset)


    # ts.find_motif()

    cnames = ["file_type","count", "expect_random","expect_nextgen", "chi_squared", "nfixed",
    "showmotif", "num", "othernum", "overlap", "ep", "ab",
    "nseqs", "v_rep_counts", "j_rep_counts"]
    motif_fn = 'tcrdist/test_files/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_cdr3_motifs_PA.log'
    x = open(motif_fn, "r").readlines()
    ts.motif_df = pd.DataFrame([l.split() for l in x], columns = cnames)

    i = 0
    row = ts.motif_df.iloc[i,:].to_dict()

    motif_list = list()
    motif_logo = list()
    for i,row in ts.motif_df.iterrows():
        StoreIOMotif_instance = ts.eval_motif(row)
        motif_list.append(StoreIOMotif_instance)
        motif_logo.append(plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600))
        if i > 1:
            break

def test_complete_example_without_motifs_step():

    import pandas as pd
    import numpy as np
    import tcrdist as td


    from tcrdist import mappers
    from tcrdist.repertoire import TCRrep
    from tcrdist.cdr3_motif import TCRMotif
    from tcrdist.subset import TCRsubset
    from tcrdist.storage import StoreIOMotif, StoreIOEntropy
    from tcrdist.plotting import plot_pwm

    tcrdist_clone_fn = 'tcrdist/test_files/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'
    tcrdist_clone_df = pd.read_csv(tcrdist_clone_fn, sep = "\t")               #1

    ind = (tcrdist_clone_df.epitope == "PA") | (tcrdist_clone_df.epitope == "F2")
    tcrdist_clone_df = tcrdist_clone_df[ind].copy()

    mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping                     #3
    tcrdist2_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df,         #4
                                                mapping = mapping)
        #1
    tr = TCRrep(cell_df = tcrdist2_df, organism = "mouse")

    #2
    tr.infer_cdrs_from_v_gene(chain = 'alpha', imgt_aligned=True)
    tr.infer_cdrs_from_v_gene(chain = 'beta',  imgt_aligned=True)

    #3
    tr.index_cols = ['clone_id', 'subject', 'epitope',
                     'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',
                     'cdr3_a_aa', 'cdr3_b_aa',
                     'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                     'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                     'cdr3_b_nucseq', 'cdr3_a_nucseq',
                     'va_countreps', 'ja_countreps',
                     'vb_countreps', 'jb_countreps',
                     'va_gene', 'vb_gene',
                     'ja_gene', 'jb_gene']

    #4
    tr.deduplicate()

    #5
    tr._tcrdist_legacy_method_alpha_beta()

    #6
    distA = tr.dist_a
    distB = tr.dist_b
    assert np.all(((distA + distB) - tr.paired_tcrdist) == 0)

    # 1
    criteria = tr.clone_df.epitope == "PA"
    clone_df_subset = tr.clone_df[criteria]

    # 2
    distA_subset = distA.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()
    distB_subset = distB.loc[clone_df_subset.clone_id, clone_df_subset.clone_id].copy()

    # 3
    ts = TCRsubset(clone_df_subset,
                 organism = "mouse",
                 epitopes = ["PA"],
                 epitope = "PA",
                 chains = ["A","B"],
                 dist_a = distA_subset,
                 dist_b = distB_subset)


    tm = TCRMotif(ts.tcr_motif_clones_df(),    #1
             organism = "mouse",
             chains = ["A"],
             epitopes = ["PA"])

    cnames = ["file_type","count", "expect_random","expect_nextgen", "chi_squared", "nfixed",
    "showmotif", "num", "othernum", "overlap", "ep", "ab",
    "nseqs", "v_rep_counts", "j_rep_counts"]

    motif_fn = 'tcrdist/test_files/mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_cdr3_motifs_PA.log'
    x = open(motif_fn, "r").readlines()
    motif_df = pd.DataFrame([l.split() for l in x], columns = cnames)
    i = 0
    row = motif_df.iloc[i,:].to_dict()

    from tcrdist.plotting import plot_pwm
    from tcrdist.storage import StoreIOMotif, StoreIOEntropy
    import pandas as pd

    # 1
    StoreIOMotif_instance = StoreIOMotif(**row)
    # 2
    StoreIOMotif_instance = ts.analyze_motif(s = StoreIOMotif_instance)
    # 3
    StoreIOMotif_instance = ts.analyze_matches(s = StoreIOMotif_instance)
    svg = plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600)
