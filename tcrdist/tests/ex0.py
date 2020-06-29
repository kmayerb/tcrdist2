import pytest

def test_ex1():
	import pandas as pd
	from tcrdist.repertoire import TCRrep
	df = pd.read_csv('dash.csv')
	df = df[df.epitope.isin(['PA'])]
	tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
	
	for chain in ['alpha','beta']:
		tr.infer_cdrs_from_v_gene(chain=chain,  imgt_aligned=True)

	tr.index_cols =['epitope','subject',
                    'cdr3_a_nucseq', 'cdr3_b_nucseq', 
                    'v_a_gene', 'j_a_gene', 
                    'v_b_gene', 'j_b_gene',
                    'cdr3_a_aa', 'pmhc_a_aa', 'cdr2_a_aa', 'cdr1_a_aa',
                    'cdr3_b_aa', 'pmhc_b_aa', 'cdr2_b_aa', 'cdr1_b_aa']

	tr.deduplicate()