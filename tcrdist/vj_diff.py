import os
import pandas as pd
from tcrsampler.sampler import TCRsampler
import numpy as np
import warnings

__all__ = ['_vj_stat', '_vj_stats']

def _vj_stats(chain = 'beta', 
              freq_type = 'vj_occur_freq', 
              df = None, 
              ts = None, 
              cols = None, 
              validate = True):
    """
    Return estimates  v-gene, j-gene, or vj-gene-pairings frequency
    specified by columns < cols > in a Pandas Dataframe < df > , given 
    a tcrsamper instance < ts >

    Make sure you provide appropriate columns with your < freq_type >. 
    For instance if freq_type = v_freq, cols must specify v-gene names column 
    only. e.g., cols = ['v_b_gene']

    Parameters
    ----------
    chain : str
        'alpha', 'beta', 'gamma', 'delta' only needed if cols are not provided.
    freq_type : str
        'vj_occur_freq', 'vj_freq', 'v_occur_freq', 'v_freq', 'j_occur_freq', 'j_freq'
    df : pd.DataFrame
        DataFrame containing v and j gene names
    ts : tcrsampler.sampler.TCRsampler
        sampler instance
    validate : bool
        If True, some basic pre-validation is performed on inputs. This is meant for 
        integration with tcrdist2. 
    cols : list
        list of column names e.g., ['v_b_gene', 'j_b_gene'] or ['v_b_gene']

    Returns
    -------
    vj_probs : list
        list of floats repressenting estimates of probability of 
        observing the specified v gene, j gene, or v,j gene pairing 
        specified in the input < df > dataframe.
    
    Examples 
    --------
    >>> import pandas as pd
    >>> import os
    >>> from tcrsampler.sampler import TCRsampler
    >>> from tcrdist.vj_diff import *
    >>> t = TCRsampler()
    >>> fn = os.path.join("tcrdist", "test_files", 'britanova_chord_blood_sample_5000.csv' )
    >>> t.ref_df = pd.read_csv(fn)
    >>> t.build_background()
    >>> dfx = pd.DataFrame({'v_b_gene':['TRBV20-1*01','TRBV20-1*01'],'j_b_gene': ['TRBJ2-1*01','TRBJ2-6*01']})
    >>> _vj_stats(chain = 'beta', df = dfx , ts = t, freq_type ='vj_occur_freq', cols = ['v_b_gene', 'j_b_gene'],  validate = True )
    [0.014802960592118424, 0.000600120024004801]
    >>> _vj_stats(chain = 'beta', df = dfx, ts = t)
    [0.014802960592118424, 0.000600120024004801]
    >>> _vj_stats(chain = 'beta', df = dfx, ts = t, cols = ['v_b_gene'], freq_type = 'v_occur_freq')
    [0.060012002400480095, 0.060012002400480095]
    >>> _vj_stats(chain = 'beta', df = dfx, ts = t, cols = ['j_b_gene'], freq_type = 'j_occur_freq')
    [0.272254450890178, 0.013802760552110422]
    """
    tcrdist2_genes =  { 'alpha' :  ['v_a_gene', 'j_a_gene'],
                        'beta'  :  ['v_b_gene', 'j_b_gene'],
                        'gamma'  : ['v_g_gene', 'j_g_gene'],
                        'delta'  : ['v_d_gene', 'j_d_gene']}
    if cols is None:
        cols = tcrdist2_genes[chain]
    
    if validate: # IF validate == True ; THEN do some basic validation of inputs.
        _vj_stats_validate_inputs(chain = chain, freq_type = freq_type, df=df, ts=ts, cols = cols)

    v_j = df[cols]
    vjprobs = list()
    for i,r in v_j.iterrows():
        try:
            if len(cols) == 2:
                assert freq_type in ['vj_freq', 'vj_occur_freq'], f"TRYING TO ACCESS {freq_type} WITH ONLY V OR J GENE PROVIDED, CHECK < cols > "
                k = tuple(r[cols].tolist())
            elif len(cols) == 1:
                assert freq_type in ['v_freq', 'v_occur_freq','j_freq', 'j_occur_freq'], f"YOU ARE TRYING TO ACCESS {freq_type} WITH BOTH V OR J GENE PROVIDED, CHECK < cols >"
                k = r[cols].tolist()[0]
            vjprob = ts.__dict__[freq_type][k]
        except KeyError:
            vjprob = 0.0
        vjprobs.append(vjprob)
    return vjprobs
     

def _vj_stats_validate_inputs(chain,freq_type, df, ts, cols):
    """
    A Pure Validation Routine for _vj_stats to ensure inputs
    are likely to work with tcrdist2 

    Parameters
    ----------
    df : pd.DataFrame
   
    ts : tcrsampler.sampler.TCRsampler
    
    Returns
    -------
    valid : bool 
        True if inputs are valid

    """
    valid = False
    tcrdist2_genes =  { 'alpha' :  ['v_a_gene', 'j_a_gene'],
                        'beta'  :  ['v_b_gene', 'j_b_gene'],
                        'gamma'  : ['v_g_gene', 'j_g_gene'],
                        'delta'  : ['v_d_gene', 'j_d_gene']}
    
    if chain not in tcrdist2_genes.keys():
        raise KeyError("chain must be 'alpha', 'beta', 'gamma', or 'delta' ")
    
    if freq_type not in ['vj_occur_freq', 'vj_freq', 'v_occur_freq', 'v_freq', 'j_occur_freq','j_freq']:
        raise KeyError("freq type must be one of 'vj_occur_freq', 'vj_freq', 'v_occur_freq', 'v_freq', 'j_occur_freq', or 'j_freq'")
    
    if freq_type in ['vj_occur_freq','vj_freq']:
        assert len(cols) == 2, f"Given {freq_type}, you must supply two columns, one each for the v and j genes"
    if freq_type in ['v_occur_freq','v_freq' ]:
        assert len(cols) == 1, f"Given {freq_type}, you must supply only one appropriate column, for the v gene"
    if freq_type in ['j_occur_freq','j_freq']:
        assert len(cols) == 1, f"Given {freq_type}, you must supply only one appropriate column, for and j gene"

    if ts is None:
        raise TypeError('_vj_stats argument < ts > must be a tcrsampler.sampler.TCRsampler instance')
    
    if not isinstance(df, pd.DataFrame): 
        raise TypeError('_vj_stats argument < df > must be a pd.DataFrame instance')

    all_tcrdist2_genes = np.concatenate(list(tcrdist2_genes.values()))
    if not np.sum([col in all_tcrdist2_genes for col in df.columns ]) >= 1:
        raise KeyError(f'_vj_stats argument < df > DataFrame must contain 1 or more columns from {all_tcrdist2_genes}')
    valid = True
    return valid

def _vj_stat(v = None, j = None,  freq_type = 'vj_occur_freq', ts = None):
    """
    Return estimate of a single  v-gene, j-gene, or vj-gene-pairings frequency
    specified < v > and <j> argumens , given a tcrsamper instance < ts >

    Parameters
    ----------
    v : str
    j : str
        e.g., 
    freq_type : str
        'vj_occur_freq', 'vj_freq', 'v_occur_freq', 'v_freq', 'j_occur_freq', 'j_freq'
    df : pd.DataFrame
        DataFrame containing v and j gene names
    ts : tcrsampler.sampler.TCRsampler
        sampler instance
    
    Example 
    -------
    >>> import pandas as pd
    >>> import os
    >>> from tcrsampler.sampler import TCRsampler
    >>> from tcrdist.vj_diff import *
    >>> t = TCRsampler()
    >>> fn = os.path.join("tcrdist", "test_files", 'britanova_chord_blood_sample_5000.csv' )
    >>> t.ref_df = pd.read_csv(fn)
    >>> t.build_background()
    >>> _vj_stat(v = 'TRBV20-1*01' , j ='TRBJ2-1*01',  ts = t, freq_type = 'vj_occur_freq')
    0.014802960592118424
    >>> _vj_stat(v = 'TRBV20-1*01' , ts = t, freq_type = 'v_occur_freq')
    0.060012002400480095
    >>> _vj_stat(j = 'TRBJ2-1*01',  ts = t, freq_type = 'j_occur_freq')
    0.272254450890178
    """
    if ts is None:
        raise ValueError("._vj_stat requires < ts > be a TCRsampler instance")
    
    if v is None and j is None:
        raise ValueError("Niether a v- nor j-gene was supplied to ._vj_stat ; atleast one must be provided")
    
    if v is None:
        tp = j
        assert freq_type in ['j_freq', 'j_occur_freq']
    elif j is None:
        tp = v
        assert freq_type in ['v_freq', 'v_occur_freq']
    else:
        tp = (v,j)
        assert freq_type in ['vj_freq', 'vj_occur_freq']

    return ts.__dict__[freq_type][tp]
