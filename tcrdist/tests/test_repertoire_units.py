from tcrdist.repertoire import _map_clone_df_to_TCRMotif_clone_df
import pytest
import pandas as pd
import numpy as np
tempSkip = pytest.mark.skip(reason="Temporarily skipping for efficiency.")

def test__map_clone_df_to_TCRMotif_clone_df():
    """
    test that _map_clone_df_to_TCRMotif_clone_df
    converts clone_df DataFrame used in tcrdist2 to the input clones_df
    DataFrame required by TCRMotif().
    """
    df = pd.DataFrame([[1,2,3,4,5,6,7,8,9,10]],
                        columns = [ 'subject','epitope',
                                    'v_a_gene', 'j_a_gene',
                                    'v_b_gene', 'j_b_gene',
                                    'cdr3_a_aa','cdr3_b_aa',
                                    'cdr2_a_aa', 'cdr2_b_aa'])
    r_df = _map_clone_df_to_TCRMotif_clone_df(df).to_dict()
    expect_df = {'subject': {0: 1}, 'epitope': {0: 2}, 'va_rep': {0: 3}, 'ja_rep': {0: 4}, 'vb_rep': {0: 5}, 'jb_rep': {0: 6}, 'cdr3a': {0: 7}, 'cdr3b': {0: 8}}
    assert r_df == expect_df

def test__map_clone_df_to_TCRMotif_clone_df_raises_KeyError():
    """
    Test that KeyError alerts that clone_df lacks correct column
    names.

    In this case: KeyError: 'clone_df must have columns names: subject'
    """
    with pytest.raises(KeyError):
        df = pd.DataFrame([[1,2,3,4,5,6,7,8,9,10]],
                            columns = [ 'X','epitope',
                                        'v_a_gene', 'j_a_gene',
                                        'v_b_gene', 'j_b_gene',
                                        'cdr3_a_aa','cdr3_b_aa',
                                        'cdr2_a_aa', 'cdr2_b_aa'])
        _map_clone_df_to_TCRMotif_clone_df(df)

def test__map_clone_df_to_TCRMotif_clone_df_KeyError_message():
    """
    Test that KeyError alerts that clone_df lacks correct column
    names.

    In this case: KeyError: 'clone_df must have columns names: subject'
    """
    with pytest.raises(KeyError) as excinfo:
        df = pd.DataFrame([[1,2,3,4,5,6,7,8,9,10]],
                            columns = [ 'X','epitope',
                                        'v_a_gene', 'j_a_gene',
                                        'v_b_gene', 'j_b_gene',
                                        'cdr3_a_aa','cdr3_b_aa',
                                        'cdr2_a_aa', 'cdr2_b_aa'])
        _map_clone_df_to_TCRMotif_clone_df(df)
    assert 'clone_df must have columns names: subject' in str(excinfo.value)

@tempSkip
def test_data_type_conversion_with_reduce_file_size():
    """
    Test that successfull conversion of np attributes from
    np.float64 to np.int16
    after calling

    """
    import pandas as pd
    import numpy as np
    import tcrdist as td

    from tcrdist import mappers
    from tcrdist.repertoire import TCRrep

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
    #print(type(tr.cdr3_a_aa_pw[1,1]))
    assert isinstance(tr.cdr3_a_aa_pw[1,1], np.int)
    assert isinstance(tr.cdr3_b_aa_pw[1,1], np.int)
    tr.reduce_file_size()
    assert isinstance(tr.cdr3_a_aa_pw[1,1], np.int16)
    assert isinstance(tr.cdr3_b_aa_pw[1,1], np.int16)


@pytest.fixture(scope="module")
def generate_tr():
    import pandas as pd
    import numpy as np
    import tcrdist as td


    from tcrdist import mappers
    from tcrdist.repertoire import TCRrep

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
    return tr
