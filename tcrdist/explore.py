# Interactive Visualization with Altair
# 2020-05-28

import altair as alt
import pandas as pd
import numpy as np
import warnings

def interactive_tsne(clone_df,
    x_var='tSNE1',
    y_var='tSNE2', 
    shape_var = 'public',
    color_var = 'log10pgen',
    color_scheme = 'magma',
    color_domain  = None, 
    mark_size = 60,
    shape_range = ['circle', 'triangle-right'],
    tooltips = ['cluster_index', 'log10pgen' ,'cluster_count', "cdr3_b_aa", "v_b_gene", "j_b_gene"]):
    """
    clone_df : pd.DataFrame,
    x_var : str
        usually 'tSNE1' must match name in clone_df
    y_var : str
        usually 'tSNE2' must match name in clone_df
    shape_var : str  
        variable for determining shape 
    color_var : str
        variable for determining color
    color_domain  : None or list
    
    
    mark_size : int
    
    shape_range : list
        list of shapes must be same length or logner than number of categories in shape_var
    tooltips :
        list of variables to include in tooltip

    """
    
    assert isinstance(clone_df, pd.DataFrame)
    assert x_var in clone_df.columns, f"{x_var} is missing from the TCRrep.clone_df instance"
    assert y_var in clone_df.columns, f"{y_var} is missing from the TCRrep.clone_df instance"
    assert color_var in clone_df.columns,  f"{color_var} sissing from the TCRrep.clone_df instance, see infer_olga_aa_cdr3_pgens"
    assert np.all([tt in clone_df.columns for tt in tooltips]), "Some of your tooltips are missing from the TCRrep.clone_df instance"
    
    shape_domain = list(clone_df[shape_var].value_counts().keys())
    assert len(shape_range) == len(shape_domain)
    #shape_domain = [{True : 'public', False: 'private'}[x] for x in shape_domain]

    print(shape_domain)

    if color_domain is None:
        some_really_small_number = 0
        if min(clone_df[color_var]) < some_really_small_number:
            warnings.warn(f"color_var has an extreme min value,  with range:{min(clone_df[color_var])}:{max(clone_df[color_var])}, you should set color_domain manually")
            dmin = 0
        else:
            dmin =  min(clone_df[color_var])
        some_really_big_number = 10000
        if max(clone_df[color_var]) > some_really_big_number:
            warnings.warn(f"color_var has an extreme max, with range: {min(clone_df[color_var])}:{max(clone_df[color_var])}, you should set color_domain manually")
            dmax = 3*np.median(clone_df[color_var])
        else:
            dmax = max(clone_df[color_var])


        color_domain = [dmin, dmax]
    
    graphic = alt.Chart(clone_df).mark_point(size=mark_size).encode(
        x=x_var,
        y=y_var,
        color= alt.Color(color_var, scale=alt.Scale(scheme=color_scheme, domain=color_domain)),
        #fill=alt.Color(color_var, scale=alt.Scale(scheme=color_scheme, domain=color_domain)),
        shape=alt.Shape('public', scale=alt.Scale(range=['triangle', 'circle'])),
        tooltip=tooltips).interactive()

    return graphic

