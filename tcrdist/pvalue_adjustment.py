import numpy as np
import pandas as pd
import statsmodels.api as sm
from functools import partial

__all__ = ['adjustnonnan']

def adjustnonnan(pvalues, method='holm'):
    """Convenient function for doing p-value adjustment.
    Accepts any matrix shape and adjusts across the entire matrix.
    Ignores nans appropriately.

    Parameters
    ----------
    pvalues : list, pd.DataFrame, pd.Series or np.ndarray
        Contains pvalues and optionally nans for adjustment.
    method : str
        An adjustment method for sm.stats.multipletests.
        Use 'holm' for Holm-Bonferroni FWER-adj and
        'fdr_bh' for Benjamini and Hochberg FDR-adj

    Returns
    -------
    adjpvalues : same as pvalues in type and shape"""

    """Turn it into a one-dimensional vector"""

    p = np.asarray(pvalues).ravel()

    """adjpvalues intialized with p to copy nans in the right places"""
    adjpvalues = np.copy(p)
    
    nanInd = np.isnan(p)
    p = p[~nanInd]
    if len(p) == 0:
        return pvalues
        
    """Drop the nans, calculate adjpvalues, copy to adjpvalues vector"""
    rej, q, alphasidak, alphabon = sm.stats.multipletests(p, alpha=0.05, method=method)
    adjpvalues[~nanInd] = q
    
    """Reshape adjpvalues"""
    if not isinstance(pvalues, list):
        adjpvalues = adjpvalues.reshape(pvalues.shape)

    """Return same type as pvalues"""
    if isinstance(pvalues, list):
        return [pv for pv in adjpvalues]
    elif isinstance(pvalues, pd.core.frame.DataFrame):
        return pd.DataFrame(adjpvalues, columns=pvalues.columns, index=pvalues.index)
    elif isinstance(pvalues, pd.core.series.Series):
        return pd.Series(adjpvalues, name=pvalues.name, index=pvalues.index)
    else:
        return adjpvalues