import os
import pytest
import pandas as pd
import numpy as np
from tcrsampler.sampler import TCRsampler
from tcrdist.vj_diff import *
import tcrdist as td
from tcrdist.repertoire import TCRrep

def _setup_tcrsampler():
    fn = os.path.join('tcrdist', 'test_files', 'britanova_chord_blood_sample_5000.csv' )
    t = TCRsampler()
    np.random.seed(11)
    t.ref_df = pd.read_csv(fn)#.sample(10000, replace=True)
    t.build_background()
    return t

t = _setup_tcrsampler()

def test_vj_stats():
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
    dfx = t.ref_df.head(10).rename(columns = {'v_reps':'v_b_gene', 'j_reps':'j_b_gene'})
    with pytest.raises(AssertionError):
        r = _vj_stats(  chain = 'beta', 
                    freq_type =ft, #'vj_occur_freq',  # <<<<-------
                    df = dfx , 
                    ts = t, 
                    cols = cols,
                    validate = True)

ft_testspace = [('TRBV20-1*01' ,'TRBJ2-1*01','vj_occur_freq',t),
                ('TRBV20-1*01' ,None,'v_occur_freq',t),
                (None, 'TRBJ2-1*01','j_occur_freq',t),
                ('TRBV20-1*01' ,'TRBJ2-1*01','vj_freq',t),
                ('TRBV20-1*01',None,'v_freq',t),
                (None,'TRBJ2-1*01','j_freq',t)]
@pytest.mark.parametrize("v,j,ft,ts", ft_testspace)
def test_vj_stat(v, j, ft, ts):
    r = _vj_stat(v=v, j=j, freq_type=ft, ts=ts)
    # print(r)
    assert isinstance(r, float)

def test_vj_surprise():
    data_fn = os.path.join('tcrdist', 'test_files', 'vdjDB_PMID28636592.tsv' )
    t_df = td.mappers.vdjdb_to_tcrdist2(pd_df=pd.read_csv(data_fn, sep='\t'))

    t_df = t_df.loc[(t_df.organism == 'HomoSapiens')]
    tr = TCRrep(cell_df=t_df, organism='human')
    tr.index_cols =['subject',
                    'cdr3_b_aa',
                    'epitope',
                    'v_b_gene',
                    'j_b_gene']
    tr.deduplicate()
       
    m1 = vj_surprise(tr.clone_df.loc[tr.clone_df.epitope == 'M1'], t, chain='beta',
                         bootstrap_samples=100)
    assert isinstance(m1, pd.Series)

"""Example plot:
data_fn = opj(_fg_data, 'tcrdist', 'datasets', 'vdjDB_PMID28636592.tsv')
t_df = td.mappers.vdjdb_to_tcrdist2(pd_df=pd.read_csv(data_fn, sep='\t'))
t_df = t_df.loc[(t_df.organism == 'HomoSapiens')]

tr = TCRrep(cell_df=t_df, organism='human')
tr.infer_cdrs_from_v_gene(chain='alpha')
tr.infer_cdrs_from_v_gene(chain='beta')
tr.index_cols =['subject',
                'cdr3_b_aa',
                'epitope',
                'v_b_gene',
                'j_b_gene']
tr.deduplicate()


    
m1 = vj_surprise(tr.clone_df.loc[tr.clone_df.epitope == 'M1'], t, chain='beta',
                     bootstrap_samples=10000)

bmlf = vj_surprise(tr.clone_df.loc[tr.clone_df.epitope == 'BMLF'], t, chain='beta',
                     bootstrap_samples=10000)

pp65 = vj_surprise(tr.clone_df.loc[tr.clone_df.epitope == 'pp65'], t, chain='beta',
                     bootstrap_samples=10000)

res = [m1, bmlf, pp65]
epitopes = ['M1', 'BMLF', 'pp65']

key = 'lr'
colors = ['C1', 'C2', 'C3']

figh = plt.figure()
axh = plt.axes((0.1, 0.1, 0.8, 0.8))
x = 0
xt = []
xlab = []
for segment in ['v', 'j', 'vj']:
    for ep, r, color in zip(epitopes, res, colors):
        if ep == 'BMLF':
            xt.append(x)
            xlab.append(segment.upper())
        if segment == 'v':
            label = ep
        else:
            label = None
        plt.plot([x, x], [r['%s_%s_ll' % (segment, key)], r['%s_%s_ul' % (segment, key)]], '-', lw=2, color=color)
        plt.plot(x, r['%s_%s' % (segment, key)], 's', color=color, label=label)
        x += 1
    x += 3
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.xticks(xt, xlab)
plt.xlabel('Segment(s)')
plt.ylabel(key.upper())
"""