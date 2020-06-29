import pytest 

def test_ex11():
	import pandas as pd
	import numpy as np
	from tcrdist.repertoire import TCRrep
	tr = TCRrep(cell_df=pd.DataFrame(), chains=['alpha','beta'], organism='mouse')
	tr.rebuild(dest_tar_name = "some_archive.tar.gz")