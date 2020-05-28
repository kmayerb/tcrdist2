# Interactive Visualization with Altair
# 2020-05-28

import altair as alt
import pandas as pd
import numpy as np
import warnings

def interactive_tsne(tr,
    x_var='tSNE1',
    y_var='tSNE2', 
    shape_var = 'public',
    color_var = 'log10pgen',
    color_scheme = 'magma',
    color_domain  = None, 
    mark_size = 60,
    shape_range = ['circle', 'triangle-right'],
    tooltips = ['cluster_index', 'log10pgen' ,'cluster_count', "cdr3_b_aa", "v_b_gene", "j_b_gene"]):
    
    assert isinstance(tr.clone_df, pd.DataFrame)
    assert x_var in tr.clone_df.columns, f"{x_var} is missing from the TCRrep.clone_df instance"
    assert y_var in tr.clone_df.columns, f"{y_var} is missing from the TCRrep.clone_df instance"
    assert color_var in tr.clone_df.columns,  f"{color_var} sissing from the TCRrep.clone_df instance, see infer_olga_aa_cdr3_pgens"
    assert np.all([tt in tr.clone_df.columns for tt in tooltips]), "Some of your tooltips are missing from the TCRrep.clone_df instance"
    
    shape_domain = list(tr.clone_df[shape_var].value_counts().keys())
    assert len(shape_range) == len(shape_domain)
    #shape_domain = [{True : 'public', False: 'private'}[x] for x in shape_domain]

    print(shape_domain)

    if color_domain is None:
        some_really_small_number = 0
        if min(tr.clone_df[color_var]) < some_really_small_number:
            warnings.warn(f"color_var has an extreme min value,  with range:{min(tr.clone_df[color_var])}:{max(tr.clone_df[color_var])}, you should set color_domain manually")
            dmin = 0
        else:
            dmin =  min(tr.clone_df[color_var])
        some_really_big_number = 10000
        if max(tr.clone_df[color_var]) > some_really_big_number:
            warnings.warn(f"color_var has an extreme max, with range: {min(tr.clone_df[color_var])}:{max(tr.clone_df[color_var])}, you should set color_domain manually")
            dmax = 3*np.median(tr.clone_df[color_var])
        else:
            dmax = max(tr.clone_df[color_var])


        color_domain = [dmin, dmax]
    
    graphic = alt.Chart(tr.clone_df).mark_point(size=mark_size).encode(
        x=x_var,
        y=y_var,
        fill=alt.Color(color_var, scale=alt.Scale(scheme=color_scheme, domain=color_domain)),
        shape=alt.Shape('public', scale=alt.Scale(range=['cross', 'circle'])),
        tooltip=tooltips).interactive()

    return graphic

