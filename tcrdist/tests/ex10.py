import pytest 

def test_ex10():
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
	
	tr.archive( dest = "some_archive", 
				dest_tar_name = "some_archive.tar.gz", 
				verbose = True, 
				use_csv = False)
