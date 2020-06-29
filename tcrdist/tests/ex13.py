import pytest 

def test_ex13():
	import pandas as pd
	import numpy as np
	from tcrdist.repertoire import TCRrep
	df = pd.read_csv('dash.csv')
	df = df[df.epitope.isin(['PA'])]
	tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
	tr.tcrdist2(processes = 1,
				metric = 'nw',
				reduce = True,
				dump = False,
				save = False)

	tr.cluster_index = tr.simple_cluster_index(pw_distances = None,
                                               method = 'complete',
                                               criterion = "distance",
                                               t = 100)
	assert len(tr.cluster_index) == tr.clone_df.shape[0]

	tr.cluster_df = tr.cluster_index_to_df(tr.cluster_index)
	