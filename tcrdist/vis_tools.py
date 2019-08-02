import matplotlib
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt



def bostock_cat_colors(color_sets = ["set3"]):

    """
    Get almost as many categorical colors as you please.

    Get more than one of the color brewer sets with ['set1' , 'set2']


    Parameters
    ----------

    sets : list
        list of color sets to return valid options are
        (set1, set2, set3, pastel1, pastel2, paired, dark, accent, category10)

    Returns
    -------

    categorical_colors : list
        list of strings (e.g. ["#e41a1c",...])

    Examples
    --------

    >>> bostock_cat_colors(['set3'])[:5]
    ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3']

    >>> bostock_cat_colors(['category10'])[:5]
    ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    Notes
    -----

    list of hex colors can be found here:
    https://observablehq.com/@d3/color-schemes


    """
    bostock = \
    {"set1"     : ["#e41a1c","#377eb8","#4daf4a","#984ea3",
                   "#ff7f00","#ffff33","#a65628","#f781bf",
                   "#999999"],
     "set2"     : ["#66c2a5","#fc8d62","#8da0cb","#e78ac3",
                   "#a6d854","#ffd92f","#e5c494","#b3b3b3"],

     "set3"     : ["#8dd3c7","#ffffb3","#bebada","#fb8072",
                   "#80b1d3","#fdb462","#b3de69","#fccde5",
                   "#d9d9d9","#bc80bd","#ccebc5","#ffed6f"],

     "pastel1"  : ["#fbb4ae","#b3cde3","#ccebc5","#decbe4",
                   "#fed9a6","#ffffcc","#e5d8bd","#fddaec",
                   "#f2f2f2"],

     "pastel2"  : ["#b3e2cd","#fdcdac","#cbd5e8","#f4cae4",
                   "#e6f5c9","#fff2ae","#f1e2cc","#cccccc"],

     "paired"   : ["#a6cee3","#1f78b4","#b2df8a","#33a02c",
                   "#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
                   "#cab2d6","#6a3d9a","#ffff99","#b15928"],

     "dark"     : ["#1b9e77","#d95f02","#7570b3","#e7298a",
                   "#66a61e","#e6ab02","#a6761d","#666666"],

     "accent"   : ["#7fc97f","#beaed4","#fdc086","#ffff99",
                   "#386cb0","#f0027f","#bf5b17","#666666"],

     "category10":["#1f77b4","#ff7f0e","#2ca02c","#d62728",
                   "#9467bd","#8c564b","#e377c2","#7f7f7f",
                   "#bcbd22","#17becf"]

    }

    l = [bostock[k] for k in color_sets]
    categorical_colors = [item for sublist in l for item in sublist]
    return categorical_colors


def cluster_viz(px,
                clone_df,
                epitopes,
                epitope_colors,
                title):
    """
    Parameters
    ----------

    px : pandas DataFrame
        DataFrame with pairwise distances row and col order
        matching clone_df DataFrame

    clone_df : pandas DataFrame
        DataFrame containing ordered list of epitope

    epitopes : list
        list of strings specifying epitope names

    epitope_colors : list
        list of string specifying hex colors

    title : string
        figure title

    Returns
    -------

    g1 : seaborn clustermap

    """
    lut = dict(zip(epitopes, epitope_colors))
    row_colors = clone_df.epitope.map(lut)

    # Cluster using seaborn
    g = sns.clustermap(data= px,
                       row_colors = row_colors,
                       col_colors = row_colors,
                       row_cluster=True,
                       col_cluster=True,
                       yticklabels=False,
                       xticklabels=False,
                      )

    # Make a Custom Legend
    for label in epitopes:
        g.ax_col_dendrogram.bar(0, 0, color=lut[label],
                                label=label, linewidth=0)

    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.legend(loc="center", ncol = 4)
    g.cax.set_position([.97, .2, .03, .45])
    g.fig.suptitle(title, fontsize=20)


def single_roc(roc, key, **kwargs):
    """
    Parameters
    ----------
    df : DataFrame
        pandas DataFrame

    Yields
    ------
    plot commands

    """
    thr = roc['thr']
    fpr = roc['fpr']
    tpr = roc['tpr']
    ind_50 = np.where(thr < .5)[0][0]

    plt.plot(roc['fpr'], roc['tpr'],
             lw=1, label='({} AUC = {})'.format(key, round(roc['roc_auc'][0],2)), **kwargs)
    #plt.plot(roc['fpr'], roc['tpr'],lw=lw)
    plt.plot([0, 1], [0, 1], color='black', lw=.5, linestyle='--')
    plt.xlim([-0.01, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(key)
    plt.legend(loc="lower right")
    plt.scatter(fpr[[ind_50]], tpr[[ind_50]], label= 'P({} > 0.5)'.format({key}))
