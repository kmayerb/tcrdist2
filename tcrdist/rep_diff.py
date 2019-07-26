import pandas as pd
import numpy as np
import itertools
import warnings

import statsmodels.api as sm
import patsy

from scipy.stats import chi2_contingency
from scipy.stats.contingency import expected_freq
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance

import fishersapi

from .pvalue_adjustment import adjustnonnan

__all__ = ['neighborhoodDiff', 'hclusterDiff']

"""TODO: Parallelize permutation test for L2-penalized logistic regression"""
'''
import parmap
import multiprocessing
pool = multiprocessing.Pool(processes=2)

def _rand_glm(n, cdf, nparams, glmParams):
    rparams = np.zeros((len(nparams), n))
    for sampi in range(n):
        randy = cdf['NBR'].sample(frac=1, replace=False).values
        rres = sm.GLM(endog=randy, **glmParams).fit_regularized(L1_wt=0, alpha=0)
        rparams[:, sampi] = rres.params
    return rparams
out = parmap.map(_rand_glm, range(nperms), cdf=cdf, nparams=2, glmParams=glmParams, pm_pool=pool, pm_chunksize=nperms // (2+1))'''


def _prep_counts(cdf, xcols, ycol='NBR', count_col=None):
    if count_col is None:
        cdf = cdf.assign(Count=1)
        count_col = 'Count'
    counts = cdf.groupby(xcols + [ycol], sort=True)[count_col].agg(np.sum).unstack(ycol).fillna(0)[[0, 1]]
    return counts

def _chi2NBR(res_df, count_cols):
    res = {'chisq':np.nan * np.zeros(res_df.shape[0]),
            'pvalue':np.nan * np.zeros(res_df.shape[0])}
    for i in range(res_df.shape[0]):
        tab = res_df[count_cols].iloc[i].values.reshape((len(count_cols) // 2, 2))
        res['chisq'][i], res['pvalue'][i], dof, expected = chi2_contingency(tab)
    return res

def _chi2_fishersNBR(counts):
    labels = []
    for rowi in counts.index.tolist():
        if type(rowi) is tuple:
            labels.append('|'.join(rowi))
        else:
            labels.append(rowi)
    res = {}
    res['chisq'], res['pvalue'], dof, expected = chi2_contingency(counts.values)
    for rowi, rowj in itertools.combinations(range(counts.shape[0]), 2):
        lab = '%s vs %s' % (labels[rowi], labels[rowj])
        # OR = ((a/b) / (c/d)) or a*d/b*c
        """It is assumed here that the number clones in the neighborhood is in col_j = 1"""
        OR, p = fishersapi.fishers_vec(counts.iloc[rowi, 1],
                                       counts.iloc[rowi, 0],
                                       counts.iloc[rowj, 1],
                                       counts.iloc[rowj, 0],
                                       alternative='two-sided')
        RR = (counts.iloc[rowi, 1] / (counts.iloc[rowi, 1] + counts.iloc[rowi, 0])) / (counts.iloc[rowj, 1] / (counts.iloc[rowj, 1] + counts.iloc[rowj, 0]))
        res.update({'RR %s' % lab: RR,
                    'OR %s' % lab: OR,
                    'pvalue %s' % lab: p})
    return res

def _fisherNBR(res_df, count_cols):
    a = res_df[count_cols[0]].values
    b = res_df[count_cols[1]].values
    c = res_df[count_cols[2]].values
    d = res_df[count_cols[3]].values

    OR, p = fishersapi.fishers_vec(a, b, c, d, alternative='two-sided')
    """It is assumed here that the number clones in the neighborhood is in col_j = 1 (i.e. b, d)"""
    RR = (b / (a + b)) / (d / (c + d))
    return {'RR':RR, 'OR':OR, 'pvalue':p}

def _glmCatNBR(df, x_cols, y_col='NBR', count_col=None, l2_alpha=0, nperms=100):
    if count_col is None:
        freq_weights = None
    else:
        freq_weights = df[count_col]

    formula = ' + '.join(['C(%s)' % c for c in x_cols])
    X = patsy.dmatrix(formula, df, return_type='dataframe')
    glmParams = dict(exog=X,
                     family=sm.families.Binomial(link=sm.families.links.logit),
                     freq_weights=freq_weights,
                     hasconst=True)
    mod = sm.GLM(endog=df[y_col].values, **glmParams)
    if l2_alpha == 0:
        res = mod.fit()
        out = {'%s_pvalue' % c:res.pvalues[c] for c in X.columns if not 'Intercept' in c}
    else:
        res = mod.fit_regularized(L1_wt=0, alpha=l2_alpha)
        rparams = np.zeros((len(res.params), nperms))
        for permi in range(nperms):
            randy = df[y_col].sample(frac=1, replace=False).values
            rres = sm.GLM(endog=randy, **glmParams).fit_regularized(L1_wt=0, alpha=l2_alpha)
            rparams[:, permi] = rres.params

        perm_values = ((np.abs(res.params[:, None]) < np.abs(rparams)).sum(axis=1) + 1) / (nperms + 1)
        out = {'%s_pvalue' % c:v for c,v in zip(X.columns, perm_values) if not 'Intercept' in c}
    
    out.update({'%s_coef' % c:res.params[c] for c in X.columns if not 'Intercept' in c})
    out.update({'%s_OR' % c:np.exp(res.params[c]) for c in X.columns if not 'Intercept' in c})
    return out

def neighborhoodDiff(clone_df, pwmat, x_cols, count_col='count', test='chi2', knn_neighbors=50, knn_radius=None, test_only=None, **kwargs):
    """Tests for association of categorical variables in x_cols with the neighborhood
    around each TCR in clone_df. The neighborhood is defined by the K closest neighbors
    using pairwise distances in pwmat, or defined by a distance radius.

    Use Fisher's exact test (test='fishers') to detect enrichment/association of the neighborhood
    with one variable.

    Tests the 2 x 2 table for each clone:

            Neighborhood
             Y   N
           ---------
         1 | a | b |
    VAR    |-------|
         0 | c | d |
           ---------

    Use the chi-squared test (test='chi2') or logistic regression (test='logistic') to detect association across multiple variables.
    Note that with sparse neighborhoods Chi-squared tests and logistic regression are unreliable. It is possible
    to pass an L2 penalty to the logistic regression using l2_alpha in kwargs, howevere this requires a permutation
    test (nperms also in kwargs) to compute a value.

    Params
    ------
    clone_df : pd.DataFrame [nclones x metadata]
        Contains metadata for each clone.
    pwmat : np.ndarray [nclones x nclones]
        Square distance matrix for defining neighborhoods
    x_cols : list
        List of columns to be tested for association with the neighborhood
    count_col : str
        Column in clone_df that specifies counts.
        Default none assumes count of 1 cell for each row.
    test : str
        Specifies Fisher's exact test ("fishers"), Chi-squared ("chi2") or
        logistic regression ("glm") for testing the association.
        Also "chi2+fishers" tests the global null using Chi2 and all pairwise
        combinations of variable values using Fisher's.
    knn_neighbors : int
        Number of neighbors to include in the neighborhood.
    knn_radius : float
        Radius for inclusion of neighbors within the neighborhood.
        Specify K or R but not both.
    test_only : None or np.ndarray
        Indices into clone_df specifying the neighborhoods for testing.
    kwargs : dict
        Arguments for the various test functions (currently only logistic
        regression, which takes l2_alpha and nperms)

    Returns
    -------
    res_df : pd.DataFrame [nclones x results]
        Results from testing the neighborhood around each clone."""
    if knn_neighbors is None and knn_radius is None:
        raise(ValueError('Must specify K or radius'))
    if not knn_neighbors is None and not knn_radius is None:
        raise(ValueError('Must specify K or radius (not both)'))

    if test == 'fishers':
        test_func = _fisherNBR
        assert len(x_cols) == 1
    elif test in ['chisq', 'chi2']:
        test_func = _chi2NBR
    elif test == 'glm':
        test_func = _glmNBR

    n = clone_df.shape[0]
    assert n == pwmat.shape[0]
    assert n == pwmat.shape[1]
    ycol = 'NBR'

    if test_only is None:
        test_only = clone_df.index
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = []
        for clonei in test_only:
            ii = np.nonzero(clone_df.index == clonei)[0][0]
            if not knn_neighbors is None:
                if knn_neighbors < 1:
                    frac = knn_neighbors
                    K = int(knn_neighbors * n)
                    print('Using K = %d (%1.0f%% of %d)' % (K, 100*frac, n))
                else:
                    K = knn_neighbors
                R = np.partition(pwmat[ii, :], knn_neighbors)[knn_neighbors]
            else:
                R = knn_radius
            y = (pwmat[ii, :] <= R).astype(float)
            K = np.sum(y)

            cdf = clone_df.assign(**{ycol:y})[[ycol, count_col] + x_cols]
            counts = _prep_counts(cdf, x_cols, ycol, count_col)

            out = {'CTS%d' % i:v for i,v in enumerate(counts.values.ravel())}

            uY = [1, 0]
            out.update({'x_col_%d' % i:v for i,v in enumerate(x_cols)})
            for i,xvals in enumerate(counts.index.tolist()):
                if type(xvals) is tuple:
                    val = '|'.join(xvals)
                else:
                    val = xvals
                out.update({'x_val_%d' % i:val,
                            'x_freq_%d' % i: counts.loc[xvals, 1] / counts.loc[xvals].sum()})
            
            out.update({'index':clonei,
                        'neighbors':list(clone_df.index[np.nonzero(y)[0]]),
                        'K_neighbors':K,
                        'R_radius':R})

            if test == 'logistic':
                glm_res = _glmCatNBR(cdf, x_cols, y_col=ycol, count_col=count_col, **kwargs)
                out.update(glm_res)
            elif test == 'chi2+fishers':
                comb_res = _chi2_fishersNBR(counts)
                out.update(comb_res)
            res.append(out)

        res_df = pd.DataFrame(res)
        if test in ['fishers', 'chi2']:
            out = test_func(res_df, count_cols=[c for c in res_df.columns if c.startswith('CTS')])
            res_df = res_df.assign(**out)

    for c in [c for c in res_df.columns if 'pvalue' in c]:
        res_df = res_df.assign(**{c.replace('pvalue', 'FWERp'):adjustnonnan(res_df[c].values, method='holm'),
                                  c.replace('pvalue', 'FDRq'):adjustnonnan(res_df[c].values, method='fdr_bh')})
    return res_df

def hclusterDiff(clone_df, pwmat, x_cols, count_col='count', test='chi2', min_n=20, method='complete', **kwargs):
    """Tests for association of categorical variables in x_cols with each cluster/node
    in a hierarchical clustering of clones with distances in pwmat.

    Use Fisher's exact test (test='fishers') to detect enrichment/association of the neighborhood
    with one variable.

    Tests the 2 x 2 table for each clone:

            Neighborhood
             Y   N
           ---------
         1 | a | b |
    VAR    |-------|
         0 | c | d |
           ---------

    Use the chi-squared test (test='chi2') or logistic regression (test='logistic') to detect association across multiple variables.
    Note that with small clusters Chi-squared tests and logistic regression are unreliable. It is possible
    to pass an L2 penalty to the logistic regression using l2_alpha in kwargs, howevere this requires a permutation
    test (nperms also in kwargs) to compute a value.

    Params
    ------
    clone_df : pd.DataFrame [nclones x metadata]
        Contains metadata for each clone.
    pwmat : np.ndarray [nclones x nclones]
        Square distance matrix for defining neighborhoods
    x_cols : list
        List of columns to be tested for association with the neighborhood
    count_col : str
        Column in clone_df that specifies counts.
        Default none assumes count of 1 cell for each row.
    test : str
        Specifies Fisher's exact test ("fishers"), Chi-squared ("chi2") or
        logistic regression ("glm") for testing the association.
        Also "chi2+fishers" tests the global null using Chi2 and all pairwise
        combinations of variable values using Fisher's.
    min_n : int
        Minimum size of a cluster for it to be tested.
    kwargs : dict
        Arguments for the various test functions (currently only logistic
        regression, which takes l2_alpha and nperms)

    Returns
    -------
    res_df : pd.DataFrame [nclusters x results]
        Results from testing each cluster.
    Z : linkage matrix [clusters, 4]
        Clustering result returned from scipy.cluster.hierarchy.linkage"""
    if test == 'fishers':
        test_func = _fisherNBR
        assert len(x_cols) == 1
    elif test in ['chisq', 'chi2']:
        test_func = _chi2NBR
    elif test == 'glm':
        test_func = _glmNBR

    n = clone_df.shape[0]
    assert n == pwmat.shape[0]
    assert n == pwmat.shape[1]
    ycol = 'NBR'

    compressedDmat = distance.squareform(pwmat)
    Z = sch.linkage(compressedDmat, method=method)

    clusters = {}
    for i, merge in enumerate(Z):
        cid = 1 + i + Z.shape[0]
        clusters[cid] = [merge[0], merge[1]]

    def _get_indices(clusters, i):
        if i <= Z.shape[0]:
            return [int(i)]
        else:
            return _get_indices(clusters, clusters[i][0]) + _get_indices(clusters, clusters[i][1])

    members = {i:_get_indices(clusters, i) for i in range(Z.shape[0] + 1, max(clusters.keys()) + 1)}

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = []
        for cid, m in members.items():
            not_m = [i for i in range(n) if not i in m]
            y = np.zeros(n)
            y[m] = 1
            
            K = np.sum(y)
            if K >= min_n and K < (n-min_n):
                R = np.max(pwmat[m, :][:, m])

                cdf = clone_df.assign(**{ycol:y})[[ycol, count_col] + x_cols]
                counts = _prep_counts(cdf, x_cols, ycol, count_col)

                out = {'CTS%d' % i:v for i,v in enumerate(counts.values.ravel())}

                uY = [1, 0]
                out.update({'x_col_%d' % i:v for i,v in enumerate(x_cols)})
                for i,xvals in enumerate(counts.index.tolist()):
                    if type(xvals) is tuple:
                        val = '|'.join(xvals)
                    else:
                        val = xvals
                    out.update({'x_val_%d' % i:val,
                                'x_freq_%d' % i: counts.loc[xvals, 1] / counts.loc[xvals].sum()})
                
                out.update({'cid':cid,
                            'members_index':list(clone_df.index[m]),
                            'members_i':m,
                            'K_neighbors':K,
                            'R_radius':R})

                if test == 'logistic':
                    glm_res = _glmCatNBR(cdf, x_cols, y_col=ycol, count_col=count_col, **kwargs)
                    out.update(glm_res)
                elif test == 'chi2+fishers':
                    comb_res = _chi2_fishersNBR(counts)
                    out.update(comb_res)
                res.append(out)

        res_df = pd.DataFrame(res)
        if test in ['fishers', 'chi2']:
            out = test_func(res_df, count_cols=[c for c in res_df.columns if c.startswith('CTS')])
            res_df = res_df.assign(**out)

    for c in [c for c in res_df.columns if 'pvalue' in c]:
        res_df = res_df.assign(**{c.replace('pvalue', 'FWERp'):adjustnonnan(res_df[c].values, method='holm'),
                                  c.replace('pvalue', 'FDRq'):adjustnonnan(res_df[c].values, method='fdr_bh')})
    return res_df