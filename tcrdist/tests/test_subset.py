import pytest
import pandas as pd
from tcrdist.subset import TCRsubset
from tcrdist.tests.my_test_subset import dist_a_subset, dist_b_subset, clone_df_subset 


#############################################
### EASY CASES : When we have paired data ###
#############################################

def test_initialization_of_TCRsubset_alpha_beta_case():
    """
    Test that we can create a TCRsubset from paired input files
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    df = clone_df_subset.iloc[0:20, :].copy()
    db = dist_b_subset.iloc[0:20, 0:20]
    da = dist_a_subset.iloc[0:20, 0:20]
    ts=TCRsubset(clone_df = df,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["A","B"],
            dist_a = da,
            dist_b = db)


def test_initialization_of_TCRsubset_alpha_beta_case_plus_motif_finding():
    """
    Test that we can create a TCRsubset from paired input files
    """
    import pytest
    import pandas as pd
    from tcrdist.subset import TCRsubset
    from tcrdist.tests.my_test_subset import dist_a_subset, dist_b_subset, clone_df_subset 
    from tcrdist.cdr3_motif import TCRMotif

    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    df = clone_df_subset.iloc[0:20, :].copy()
    db = dist_b_subset.iloc[0:20, 0:20]
    da = dist_a_subset.iloc[0:20, 0:20]
    ts=TCRsubset(clone_df = df,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["A","B"],
            dist_a = da,
            dist_b = db)
    motif_df = ts.find_motif()
    assert isinstance(motif_df, pd.DataFrame)
    assert isinstance(ts.motif_df, pd.DataFrame)

  

###############################################
### HARDER CASES : When we have single data ###
###############################################

def test_initialization_of_TCRsubset_beta_case():
    """
    Test that we can create a TCRsubset from just beta chain files.
    For this to work, fake alpha sequences will have to be created:
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    df = clone_df_subset[['clone_id', 'subject', 'epitope', 
                        'v_b_gene', 'j_b_gene', 
                        'cdr3_b_aa', 'cdr1_b_aa', 'cdr2_b_aa', 
                        'pmhc_b_aa', 'cdr3_b_nucseq', 'count',
                        'vb_countreps', 'jb_countreps','vb_gene', 'jb_gene']]
    
    ts = TCRsubset(  clone_df = df,            
                organism = "mouse",
                epitopes = ["PA"] ,
                epitope = "PA",
                chains = ["beta"],
                dist_b = dist_b_subset)
    
    assert(isinstance(ts.clone_df.va_gene, pd.Series))
    assert(ts.clone_df.va_gene.iloc[0] == "TRAV10*01")
    


def test_initialization_of_TCRsubset_alpha_case():
    """
    Test that we can create a TCRsubset from paired input files
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    df = clone_df_subset[['clone_id', 'subject', 'epitope', 
                        'v_a_gene', 'j_a_gene', 
                        'cdr3_a_aa', 'cdr1_a_aa', 'cdr2_a_aa', 
                        'pmhc_a_aa', 'cdr3_a_nucseq', 'count',
                        'va_countreps', 'ja_countreps','va_gene', 'ja_gene']]
    ts = TCRsubset(  clone_df = df,            
                organism = "mouse",
                epitopes = ["PA"] ,
                epitope = "PA",
                chains = ["alpha"],
                dist_a = dist_a_subset)

    assert(isinstance(ts.clone_df.vb_gene, pd.Series))
    assert(ts.clone_df.vb_gene.iloc[0] == "TRBV1*01")


def test_initialization_of_TCRsubset_alpha_only_find_motifs():
    """
    Test that we can create a TCRsubset from paired input files
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)

    df = clone_df_subset[['clone_id', 'subject', 'epitope', 
                        'v_a_gene', 'j_a_gene', 
                        'cdr3_a_aa', 'cdr1_a_aa', 'cdr2_a_aa', 
                        'pmhc_a_aa', 'cdr3_a_nucseq', 'count',
                        'va_countreps', 'ja_countreps','va_gene', 'ja_gene']]
    df = df.iloc[0:20, :]
    
    da = dist_a_subset.iloc[0:20, 0:20]
    print(dist_a_subset)
    ts=TCRsubset(clone_df = df,            
                organism = "mouse",
                epitopes = ["PA"] ,
                epitope = "PA",
                chains = ["alpha"],
                dist_a = da)
    motif_df = ts.find_motif()
    assert isinstance(motif_df, pd.DataFrame)
    assert isinstance(ts.motif_df, pd.DataFrame)




def test_initialization_of_TCRsubset_beta_only_find_motifs():
    """
    Test that we can create a TCRsubset from paired input files
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    df = clone_df_subset[['clone_id', 'subject', 'epitope', 
                        'v_b_gene', 'j_b_gene', 
                        'cdr3_b_aa', 'cdr1_b_aa', 'cdr2_b_aa', 
                        'pmhc_b_aa', 'cdr3_b_nucseq', 'count',
                        'vb_countreps', 'jb_countreps','vb_gene', 'jb_gene']]
    df = df.iloc[0:20, :]
    
    db = dist_b_subset.iloc[0:20, 0:20]

    ts=TCRsubset(clone_df = df,            
                organism = "mouse",
                epitopes = ["PA"] ,
                epitope = "PA",
                chains = ["beta"],
                dist_b = db)

    motif_df = ts.find_motif()
    assert isinstance(motif_df, pd.DataFrame)
    assert isinstance(ts.motif_df, pd.DataFrame)



def test_my_subset_files_exist():
    """
    See my_test_subset for test files written directly in dictionary syntax
    Test that they exist
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)



@pytest.mark.skip()
def test_initialization_of_TCRsubset_alpha_beta_case():
    """
    Test that we can create a TCRsubset from paired input files
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    TCRsubset(clone_df = clone_df_subset,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["A","B"],
            dist_a = dist_a_subset,
            dist_b = dist_b_subset)





def test_initialization_of_TCRsubset_alpha_beta_case_chain_names():
    """
    Test that we can create a TCRsubset from paired input files with 
    alternate chain names
    """
    assert isinstance(dist_a_subset, pd.DataFrame)
    assert isinstance(dist_b_subset, pd.DataFrame)
    assert isinstance(clone_df_subset, pd.DataFrame)
    TCRsubset(clone_df = clone_df_subset,            
            organism = "mouse",
            epitopes = ["PA"] ,
            epitope = "PA",
            chains = ["alpha","beta"],
            dist_a = dist_a_subset,
            dist_b = dist_b_subset)


