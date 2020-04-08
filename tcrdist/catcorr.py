import pandas as pd
import numpy as np
import itertools

# from scipy.stats import chi2_contingency
# from scipy.stats.contingency import expected_freq

import fishersapi

from .pvalue_adjustment import adjustnonnan

__all__ = ['catcorr']


def catcorr(df, x_cols=None, y_cols=None, count_col=None, min_n=10):
    """Test for associations between categorical variables (columns)
    in df by testing pairs of values within pairs of columns using
    Fisher's exact test. This is not the best way to model associations
    with multinomial distributed variables, but it will work as an initial screen.

    Parameters
    ----------
    df : pd.DataFrame
        Data contains columns of categorical variables.
    min_n : int
        Minimum number of counts required for testing.
    counts_col : str
        Column name indicating counts for each row. If None
        will use one count per row.

    Returns
    -------
    res : pd.DataFrame
        Results, one row per test, with multiplicity adjustment"""
    if count_col is None:
        count_col = ''
        counts = np.ones(df.shape[0])
    else:
        counts = df[count_col].values

    if x_cols is None and y_cols is None:
        col_pairs = [p for p in itertools.combinations([c for c in df.columns if not c == count_col], 2)]
    elif x_cols is None:
        col_pairs = [p for p in itertools.combinations(y_cols, 2)]
    elif y_cols is None:
        col_pairs = [p for p in itertools.combinations(x_cols, 2)]
    else:
        col_pairs = [p for p in itertools.product(x_cols, y_cols)]

    res = []
    for col1, col2 in col_pairs:
        for val1, val2 in itertools.product(df[col1].unique(), df[col2].unique()):
            aind = (df[col1] == val1) & (df[col2] == val2)
            bind = (df[col1] == val1) & (df[col2] != val2)
            cind = (df[col1] != val1) & (df[col2] == val2)
            dind = (df[col1] != val1) & (df[col2] != val2)
            n = counts.sum()
            w = np.sum(aind.astype(int).values * counts)
            if w > min_n:
                #OR, pvalue = test_edge(df, (col1, val1), (col2, val2))
                tmp = {'xcol':col1,
                       'xval':val1,
                       'ycol':col2,
                       'yval':val2,
                       'X+Y+':w,
                       'X+Y-':np.sum(bind.astype(int).values * counts),
                       'X-Y+':np.sum(cind.astype(int).values * counts),
                       'X-Y-':np.sum(dind.astype(int).values * counts)}
                tmp.update({'X_marg':(tmp['X+Y+'] + tmp['X+Y-']) / n,
                            'Y_marg':(tmp['X+Y+'] + tmp['X-Y+']) / n,
                            'X|Y+':tmp['X+Y+'] / (tmp['X+Y+'] + tmp['X-Y+']),
                            'X|Y-':tmp['X+Y-'] / (tmp['X+Y-'] + tmp['X-Y-']),
                            'Y|X+':tmp['X+Y+'] / (tmp['X+Y+'] + tmp['X+Y-']),
                            'Y|X-':tmp['X-Y+'] / (tmp['X-Y+'] + tmp['X-Y-'])})
                res.append(tmp)

    res = pd.DataFrame(res)
    for k in ['X+Y+', 'X+Y-', 'X-Y+', 'X-Y-']:
        res = res.assign(**{k:res[k].astype(int)})

    OR, p = fishersapi.fishers_vec(res['X+Y+'],
                                   res['X+Y-'],
                                   res['X-Y+'],
                                   res['X-Y-'],
                                   alternative='two-sided')
    RR = ((res['X+Y+'] + 0.5) / (res['X+Y+'] + res['X-Y+'] + 1)) / ((res['X+Y-'] + 0.5) / (res['X+Y-'] + res['X-Y-'] + 1))
    res = res.assign(RR=RR,
                     OR=OR,
                     pvalue=p,
                     FWERp=adjustnonnan(p, method='holm'),
                     FDRq=adjustnonnan(p, method='fdr_bh'))
    res = res.sort_values(by='pvalue', ascending=True)
    cols = ['xcol', 'xval', 'ycol', 'yval',
            'RR', 'OR', 'pvalue',
            'FWERp', 'FDRq',
            'X+Y+', 'X+Y-', 'X-Y+', 'X-Y-',
            'X_marg', 'Y_marg', 'X|Y+', 'X|Y-',
            'Y|X+', 'Y|X-']
    return res[cols]