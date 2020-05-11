"""
This test page is for unit tests associated with TCRrep operations. 
As opposed to test_repertoire.py these should be exceptionally short unit tests. 


"""
from tcrdist.repertoire import TCRrep
import pytest
import pandas as pd
import numpy as np

def test_TCRrep_deduplicate_gives_warning_above_complete_values():
    """
    Warn user if there is missing value in an index column and cell count will not match clone count
    """
    df= pd.DataFrame({"cdr3_a_aa" : ["A","B","C"], 
                        "v_a_gene"  : ["TRAV1*01","TRAV1*01",None],
                        "count"     : [1,1,1]}) 
    tr = TCRrep(cell_df = df, organism = "human", chains = ["alpha"])
    tr.index_cols = ['cdr3_a_aa', 'v_a_gene'] 
    with pytest.warns(None) as record:
        tr.deduplicate()
    assert len(record) == 4
    assert str(record[0].message) == "Not all cells/sequences could be grouped into clones.\n"
    
def test_TCRrep_show_incomplete_entries():
    """
    Test that show TCRrep method returns rows with some missing attribute
    """
    df= pd.DataFrame({"cdr3_a_aa" : ["A","B","C"], 
                      "v_a_gene"  : ["TRAV1*01","TRAV1*01",None],
                      "count"     : [1,1,1]}) 
    tr = TCRrep(cell_df = df, organism = "human", chains = ["alpha"])
    tr.index_cols = ['cdr3_a_aa', 'v_a_gene'] 
    dfi = tr.show_incomplete()
    assert isinstance(dfi, pd.DataFrame)
    assert dfi.to_dict() ==  {'cdr3_a_aa': {2: 'C'}, 'v_a_gene': {2: None}}  

def test_TCRrep_show_incomplete_entries_when_there_are_none():
    """
    Test that show TCRrep method returns rows with some missing attribute
    """
    df= pd.DataFrame({"cdr3_a_aa" : ["A","B","C"], 
                      "v_a_gene"  : ["TRAV1*01","TRAV1*01","TRAV1*01"],
                      "count"     : [1,1,1]}) 
    tr = TCRrep(cell_df = df, organism = "human", chains = ["alpha"])
    tr.index_cols = ['cdr3_a_aa', 'v_a_gene'] 
    dfi = tr.show_incomplete()
    assert isinstance(dfi, pd.DataFrame)
    assert dfi.shape[0] == 0
    assert dfi.shape[1] == 2
    #assert dfi.to_dict() ==  {'cdr3_a_aa': {2: 'C'}, 'v_a_gene': {2: None}}   
    
