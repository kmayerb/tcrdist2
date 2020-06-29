import pytest 

def test_ex11():
	import pandas as pd
	from tcrsampler.sampler import TCRsampler
	fn = 'britanova_chord_blood.csv' 
	t = TCRsampler()
	t.ref_df = pd.read_csv(fn)
	t.build_background()
	t.v_freq
	t.j_freq
	t.vj_freq
	t.sample_background(v ='TRBV10-1*01', j ='TRBJ1-1*01',n=3, depth = 1, seed =1, use_frequency= True )
