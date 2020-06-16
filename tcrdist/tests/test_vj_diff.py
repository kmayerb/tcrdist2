import os
import pytest
import pandas as pd
from tcrsampler.sampler import TCRsampler
from tcrdist.vj_diff import *

def test_vj_stats():
    fn = os.path.join("tcrdist", "test_files", 'britanova_chord_blood_sample_5000.csv' )
    t = TCRsampler()
    t.ref_df = pd.read_csv(fn).sample(10000, replace = True)
    t.build_background()
   
    with pytest.raises(KeyError):
        """ column naems in < df > are not expecte """
        dfx = t.ref_df.head(10)
        _vj_stats(df = dfx , ts = t, validate = True)
    
    with pytest.raises(TypeError):
        """ Missing < df > """
        _vj_stats(ts = t, validate = True)

    with pytest.raises(TypeError):
        """ < df > is not type pd.DataFrame """
        _vj_stats(df = 1 , ts = t, validate = True)
    
    with pytest.raises(TypeError):
        """ Missing < ts > """
        dfx = t.ref_df.head(10).rename(columns = {'v_reps':'v_b_gene', 'j_reps':'j_b_gene'})
        _vj_stats(df = dfx, validate = True)
    
    dfx = t.ref_df.head(10).rename(columns = {'v_reps':'v_b_gene', 'j_reps':'j_b_gene'})
    _vj_stats(df = dfx , ts = t, validate = True)



ft_testspace = [('vj_occur_freq', ['v_b_gene', 'j_b_gene']),
                ('v_occur_freq', ['v_b_gene']),
                ('j_occur_freq', ['j_b_gene']),
                ('vj_freq', ['v_b_gene', 'j_b_gene']),
                ('v_freq', ['v_b_gene']),
                ('j_freq',['j_b_gene'])]
@pytest.mark.parametrize("ft,cols", ft_testspace)
def test_vj_stats_on_multiple_metrics(ft, cols):
    fn = os.path.join("tcrdist", "test_files", 'britanova_chord_blood_sample_5000.csv' )
    t = TCRsampler()
    t.ref_df = pd.read_csv(fn).sample(10000, replace = True)
    t.build_background()
    dfx = t.ref_df.head(10).rename(columns = {'v_reps':'v_b_gene', 'j_reps':'j_b_gene'})
    r = _vj_stats(  chain = 'beta', 
                    freq_type =ft, #'vj_occur_freq',  # <<<<-------
                    df = dfx , 
                    ts = t, 
                    cols = cols,
                    validate = True)
    print(r)
    assert isinstance(r, list)
    assert isinstance(r[0], float)

ft_testspace = [ ('j_freq',['v_b_gene','j_b_gene']),
                 ('v_freq',['v_b_gene','j_b_gene']),
                 ('vj_freq',['j_b_gene'])]
@pytest.mark.parametrize("ft,cols", ft_testspace)
def test_vj_stats_on_multiple_metrics_asserterror(ft, cols):
    fn = os.path.join("tcrdist", "test_files", 'britanova_chord_blood_sample_5000.csv' )
    t = TCRsampler()
    t.ref_df = pd.read_csv(fn).sample(10000, replace = True)
    t.build_background()
    dfx = t.ref_df.head(10).rename(columns = {'v_reps':'v_b_gene', 'j_reps':'j_b_gene'})
    with pytest.raises(AssertionError):
        r = _vj_stats(  chain = 'beta', 
                    freq_type =ft, #'vj_occur_freq',  # <<<<-------
                    df = dfx , 
                    ts = t, 
                    cols = cols,
                    validate = True)


fn = os.path.join("tcrdist", "test_files", 'britanova_chord_blood_sample_5000.csv' )
t = TCRsampler()
t.ref_df = pd.read_csv(fn).sample(10000, replace = True)
t.build_background()
ft_testspace = [('TRBV20-1*01' ,'TRBJ2-1*01','vj_occur_freq',t),
                ('TRBV20-1*01' ,None,'v_occur_freq',t),
                (None, 'TRBJ2-1*01','j_occur_freq',t),
                ('TRBV20-1*01' ,'TRBJ2-1*01','vj_freq',t),
                ('TRBV20-1*01',None,'v_freq',t),
                (None,'TRBJ2-1*01','j_freq',t)]
@pytest.mark.parametrize("v,j,ft,ts", ft_testspace)
def test_vj_stat(v,j,ft,ts):
    r = _vj_stat(v = v ,j = j, freq_type= ft, ts = ts)
    print(r)
    assert isinstance(r, float)
