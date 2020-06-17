import os
import pandas as pd
import numpy as np
import warnings

from scipy import stats
from scipy import special
from scipy.stats import multinomial

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    """Seems to be sklearn version dependent: this is consistent with the v0.23 API"""
    from sklearn.metrics import adjusted_mutual_info_score
    from sklearn.metrics.cluster.supervised import expected_mutual_information
    from sklearn.metrics.cluster import mutual_info_score
    from sklearn.metrics.cluster.supervised import _generalized_average

from tcrsampler.sampler import TCRsampler

__all__ = ['vj_surprise',
           '_vj_stat',
           '_vj_stats']

"""TODO
 - Example notebook reproducing the original Dash et al. heatmaps?
 - Consider idea of releasing numba-based bootstrap code and using numba for bootstrapping
   Depends on how useful these CIs turn out to be"""

def _js(a, b):
    """Jensen-Shannon distance
    symetrized and smoothed KL-divergence"""
    m = (a + b)
    m /= 2.
    m = np.where(m, m, 1.)
    return 0.5 * np.sum(special.xlogy(a, a/m) + special.xlogy(b, b/m))

def _entropy(a):
    """Shannon entropy"""
    return -np.sum(special.xlogy(a, a))

def _kl(a, b):
    """Kullback-Liebler divergence"""
    return np.sum(special.xlogy(a, a / b))

def _mi(a, b, ab):
    """Mutual information"""
    return np.sum(special.xlogy(ab, ab / np.dot(a[:, None], b[None, :])))

def _ami(ab_cts, average_method='arithmetic'):
    """Adjusted mutual information between two discrete categorical random variables
    based on counts observed and provided in ab_cts.

    Code adapted directly from scikit learn AMI to
    accomodate having counts/contingency table instead of rows/instances:
    https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html

    Parameters
    ----------
    ab_cts : np.ndarray [len(a_classes) x len(b_classes)
        Counts for each combination of classes in random variables a and b
        organized in a rectangular array.
    average_method : str
        See sklearn documentation for details

    Returns
    -------
    ami : float
        Adjusted mutual information score for variables a and b"""
    a_freq = np.sum(ab_cts, axis=1)
    a_freq = a_freq / np.sum(a_freq)
    b_freq = np.sum(ab_cts, axis=0)
    b_freq = b_freq / np.sum(b_freq)
    n_samples = np.sum(ab_cts)

    """ Calculate the MI for the two clusterings
    contingency is a joint count distribution [a_classes x b_classes]"""
    mi = mutual_info_score(None, None, contingency=ab_cts)

    """Calculate the expected value for the mutual information"""
    emi = expected_mutual_information(ab_cts, n_samples)
    """Calculate entropy"""
    h_true, h_pred = _entropy(a_freq), _entropy(b_freq)
    normalizer = _generalized_average(h_true, h_pred, average_method)
    denominator = normalizer - emi

    if denominator < 0:
        denominator = min(denominator, -np.finfo('float64').eps)
    else:
        denominator = max(denominator, np.finfo('float64').eps)
    ami = (mi - emi) / denominator
    return ami

def vj_surprise(member_df, tsamp, chain='beta', count_col=None,
                          v_col=None, j_col=None,
                          use_refcounts=False, bootstrap_samples=0, alpha=0.05):
    """Computes statistics about the frequency of V and J gene usage in a 
    cluster of TCR clones, relative to a reference dataset. Operates on a single
    TCR chain and can use counts of clonal expansion in the computation or not.

    These statistics are useful when attempting to identify interesting and
    antigen-specific clusters. We generally expect that antigen-specific
    clusters may be enriched for specific V, J and V-J genes that are
    required for antigen recognition. Together with similarly derived statistics
    about the CDR3, this can help rank clusters for biological validation.

    Statistics:
        Jensen-Shannon divergence (jsd): entropy-normalized JS distance between the
                TCRs in the cluster and the specified reference dataset. This
                was used in the Dash et al. (2017) manuscript. A greater
                divergence from the reference may indicate a "surprising" cluster.

        Kullback-Liebler divergence (kl): Non-symetric measure of the difference
                between two discrete frequency distributions, in this case
                the V/J frequencies observed in the cluster and those in the
                reference. A greater divergence from the reference may indicate
                a "surprising" cluster.

        Multinomial log-likelihood (loglik): Probability of observing the V/J usage in
                the cluster given the V/J usage in the reference (null hypothesis).
                The probability is related to the JS and KL divergences, but
                also takes into account the number of sequences in the cluster,
                with more TCRs divergent from the reference leading to a
                lower probability.

        Log-likelihood ratio (lr): the log of the ratio of the probability of
                observing the V/J usage under the null hypothesis (i.e. the reference)
                vs. the probability under an empirical multinomial fit to the
                cluster data. A log-LR close to 0 signifies an uninteresting
                V/J usage, while a large negative value indicates a difference
                from the reference.

    Note that for each statistic an estimate is given for the V, J and joint V-J
    gene usage in the cluster (e.g. v_jsd, j_jsd, vj_jsd. Optionally, 1 - alpha-level
    confidence intervals are provided by bootstrap resampling from member_df with replacement.

    A pseudocount is used for the reference distribution to accomodate V, J or V-J
    genes that are not observed. The pseudocount is half the minimum observed frequency.

    Also, the adjusted mutal information (AMI) between V and J gene frequencies,
    as originally reported in Dash et al. (2017), is computed as a way of
    quantifying the degree of preference and correlation of V-J genes within the cluster.

    Parameters
    ----------
    member_df : pd.DataFrame with columns v_col, j_col and count_col
        Subset of a TCR dataset, possibly a cluster for studying V/J gene usage
    tsamp : TCRSampler object
        An object that represents a reference dataset for comparison.
        See tcrsampler package documentation for details.
    chain : str
        One of alpha, beta, gamma or delta to specify which chain should be analyzed.
        If specified it will try to use default column names (e.g. v_b_gene) in
        member_df and the TCRSampler. If None, it will use v_col and j_col.
    count_col : str
        Optional, for specifying a column in member_df that contains expansion
        counts that should be used in computing gene frequencies.
    v_col, j_col : str
        Column names in member_df and TRCSampler for V and J genes.
        If None, these will default, based on the chain, to e.g. v_b_gene.
    use_refcounts : bool
        Indicator of whether the counts of expanded clones in the reference dataset
        should be used for computing reference frequencies.
    bootstrap_samples : int
        For values > 0, the number of non-parametric bootstrap samples that should
        be drawn for estimating confidence intervals on each statistic. At least 10K
        samples are recommended for good estimation.
    alpha : float [0, 1]
        If bootstrap confidence intervals are requested, the 1 - alpha-level interval
        will be provided as "_ll" and "_ul" columns for lower and upper-levels.

    Returns
    -------
    out : pd.Series
        Each element of the output contains a float value for one of the statistics
        described above.
    """
    if v_col is None or j_col is None:
        tcrdist2_genes =  { 'alpha' :  ['v_a_gene', 'j_a_gene'],
                            'beta'  :  ['v_b_gene', 'j_b_gene'],
                            'gamma'  : ['v_g_gene', 'j_g_gene'],
                            'delta'  : ['v_d_gene', 'j_d_gene']}
        v_col, j_col = tcrdist2_genes[chain]

    if count_col is None:
        v_cts = member_df.groupby(v_col).count().iloc[:, 0]
        j_cts = member_df.groupby(j_col).count().iloc[:, 0]
        vj_cts_obs = member_df.groupby([v_col, j_col]).count().iloc[:, 0]
    else:
        v_cts = member_df.groupby(v_col)[count_col].count()
        j_cts = member_df.groupby(j_col)[count_col].count()
        vj_cts_obs = member_df.groupby([v_col, j_col])[count_col].count()

    v_genes = v_cts.index.tolist()
    j_genes = j_cts.index.tolist()

    vj_ct_mat = np.zeros((v_cts.shape[0], j_cts.shape[0]))
    for (vgene, jgene), ct in vj_cts_obs.iteritems():
        vj_ct_mat[v_genes.index(vgene), j_genes.index(jgene)] = ct
    vj_freq = vj_ct_mat / np.sum(vj_ct_mat)
    vj_cts = vj_ct_mat.ravel()

    if use_refcounts:
        suffix = '_freq'
    else:
        suffix = '_occur_freq'

    v_ref = _vj_stats(chain=chain, 
                      freq_type='v' + suffix, 
                      df=v_cts.reset_index(), 
                      ts=tsamp,
                      cols=[v_col])
    
    ref_pseudocount = np.min(list(tsamp.v_freq.values())) / 2
    v_ref = (v_ref + ref_pseudocount) / np.sum((v_ref + ref_pseudocount))

    j_ref = _vj_stats(chain=chain, 
                      freq_type='j' + suffix, 
                      df=j_cts.reset_index(), 
                      ts=tsamp,
                      cols=[j_col])

    ref_pseudocount = np.min(list(tsamp.j_freq.values())) / 2
    j_ref = (j_ref + ref_pseudocount) / np.sum((j_ref + ref_pseudocount))

    vj_tmp = _vj_stats(chain=chain, 
                      freq_type='vj' + suffix, 
                      df=vj_cts_obs.reset_index(), 
                      ts=tsamp,
                      cols=[v_col, j_col])

    ref_pseudocount = np.min(list(tsamp.vj_freq.values())) / 2
    vj_ref = ref_pseudocount * np.ones((v_cts.shape[0], j_cts.shape[0]))
    for (vgene, jgene), ct in zip(vj_cts_obs.index, vj_tmp):
        vj_ref[v_genes.index(vgene), j_genes.index(jgene)] += ct
    vj_ref = vj_ref / np.sum(vj_ref)
    
    def _compute_all(vj_ct_mat):
        j_cts = vj_ct_mat.sum(axis=0)
        v_cts = vj_ct_mat.sum(axis=1)

        v_freq = v_cts / np.sum(v_cts)
        j_freq = j_cts / np.sum(j_cts)
        vj_freq = vj_ct_mat / np.sum(vj_ct_mat)

        """Legacy Dash et al."""
        v_leg = _js(v_freq, v_ref) / np.mean([_entropy(v_freq), _entropy(v_ref)])
        j_leg = _js(j_freq, j_ref) / np.mean([_entropy(j_freq), _entropy(j_ref)])
        vj_leg = _js(vj_freq.ravel(), vj_ref.ravel()) / np.mean([_entropy(vj_freq.ravel()), _entropy(vj_ref.ravel())])

        vj_ami = _ami(vj_ct_mat)
        """Below is equivalent if there are is no counts column"""
        # vj_ami = adjusted_mutual_info_score(member_df[v_col], member_df[j_col])

        """KL divergence"""
        v_kl = _kl(v_freq, v_ref)
        j_kl = _kl(j_freq, j_ref)
        vj_kl = _kl(vj_freq.ravel(), vj_ref.ravel())

        """Multinomial likelihood under the null/reference"""
        v_loglik = multinomial.logpmf(v_cts, np.sum(v_cts), v_ref)
        j_loglik = multinomial.logpmf(j_cts, np.sum(j_cts), j_ref)
        """This treats each combination of VJ as a single class"""
        vj_loglik = multinomial.logpmf(vj_ct_mat.ravel(), np.sum(vj_ct_mat), vj_ref.ravel())

        """Empirical multinomial likelihood"""
        v_emp_loglik = multinomial.logpmf(v_cts, np.sum(v_cts), v_freq)
        j_emp_loglik = multinomial.logpmf(j_cts, np.sum(j_cts), j_freq)
        """This treats each combination of VJ as a single class"""
        vj_emp_loglik = multinomial.logpmf(vj_ct_mat.ravel(), np.sum(vj_ct_mat), vj_freq.ravel())

        """Log likelihood-ratio (LR): will be near zero if cluster is distributed similar to the reference"""
        v_lr = v_loglik - v_emp_loglik
        j_lr = j_loglik - j_emp_loglik
        vj_lr = vj_loglik - vj_emp_loglik

        out = dict(v_jsd=v_leg,
                    j_jsd=j_leg,
                    vj_jsd=vj_leg,
                    vj_ami=vj_ami,
                    v_kld=v_kl,
                    j_kld=j_kl,
                    vj_kld=vj_kl,
                    v_loglik=v_loglik,
                    j_loglik=j_loglik,
                    vj_loglik=vj_loglik,
                    v_lr=v_lr,
                    j_lr=j_lr,
                    vj_lr=vj_lr)
        return out
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
    
        out = _compute_all(vj_ct_mat)
    
        if bootstrap_samples:
            samps = {}
            for k in out.keys():  
                samps[k] = np.nan * np.zeros(bootstrap_samples)
            vj_ctsr = multinomial(np.sum(vj_cts), vj_freq.ravel()).rvs(size=bootstrap_samples)
            for sampi in range(bootstrap_samples):
                outr = _compute_all(vj_ctsr[sampi, :].reshape(vj_ct_mat.shape))
                for k in outr.keys():  
                    samps[k][sampi] = outr[k]
            nvals = np.round((bootstrap_samples - 1) * np.array([alpha, 1-alpha])).astype(int)
            for k in list(out.keys()):  
                samps[k].sort()
                out[k + '_ll'] = samps[k][nvals[0]]
                out[k + '_ul'] = samps[k][nvals[1]]

    return pd.Series(out).sort_index()

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
