import pytest 

def test_ex7():
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
	tr.custom_dmat(cdr = 'cdr3_b_aa', metric =Levenshtein.distance, processes = 1)
	tr.add_custom_dmat(cdr = 'cdr3_b_aa', metric =Levenshtein.distance, processes = 1)