import pytest 

def test_ex8():
	import pandas as pd
	import numpy as np
	from tcrdist.repertoire import TCRrep
	df = pd.read_csv('dash.csv')
	df = df[df.epitope.isin(['PA'])]
	tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
	tr.tcrdist2(processes = 1,
				metric = 'hamming',
				reduce = True,
				dump = False,
				save = False)

	import Levenshtein
	for cdr,w in {'cdr3_b_aa':3,'pmhc_b_aa':1,'cdr2_b_aa':1,'cdr1_b_aa':1}.items():
		tr.add_custom_dmat(cdr = cdr, metric =Levenshtein.distance, processes = 1)

	tr.pw_levenshtein_tcrdist_beta = tr.custom_cdr3_b_aa_pw + \
	                                 tr.custom_pmhc_b_aa_pw + \
	                                 tr.custom_cdr2_b_aa_pw + \
	                                 tr.custom_cdr1_b_aa_pw 
