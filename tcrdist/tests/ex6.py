import pytest

def test_ex6():
	import pwseqdist as pw 
	from scipy.spatial.distance import squareform
	import Levenshtein
	dvec = pw.apply_pairwise_sq(seqs = ['homer', 'home', 'rome'], 
                                metric = Levenshtein.distance, 
                                ncpus = 1)
	dmat = squareform(dvec)

