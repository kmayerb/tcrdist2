import pytest

def test_TCRrep_import():
    from tcrdist.repertoire import TCRrep

def test_TCRSubet_import():
    from tcrdist.subset import TCRsubset

def test_TCRMotif_import():
    from tcrdist.cdr3_motif import TCRMotif

def test_TCRChain_import():
    from tcrdist.objects import TCRChain

def test_TCRClone_import():
    from tcrdist.objects import TCRClone

def test_StoreIO_import():
    from tcrdist.storage import StoreIOMotif, StoreIOEntropy

def test_plot_pwm_import():
    from tcrdist.plotting import plot_pwm
