import pytest

def test_ex3():
	import pandas as pd
	import numpy as np
	from tcrdist.repertoire import TCRrep
	df = pd.read_csv('dash.csv')
	df = df[df.epitope.isin(['PA'])]
	tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
	tr.tcrdist2(processes = 2,
				metric = 'hamming',
				reduce = True,
				dump = False,
				save = False)

	assert np.all( tr.pw_tcrdist == tr.pw_alpha + tr.pw_beta)