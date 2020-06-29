import pytest

def test_ex4():
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

	assert np.all( tr.pw_tcrdist == tr.pw_alpha + tr.pw_beta)
	assert np.all( tr.pw_beta    == tr.cdr3_b_aa_pw + tr.pmhc_b_aa_pw + tr.cdr2_b_aa_pw + tr.cdr1_b_aa_pw)
	assert np.all( tr.pw_alpha   == tr.cdr3_a_aa_pw + tr.pmhc_a_aa_pw + tr.cdr2_a_aa_pw + tr.cdr1_a_aa_pw)
