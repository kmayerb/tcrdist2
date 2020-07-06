from collections import namedtuple
import os
import palmotif
import pandas as pd
import pwseqdist as pw
from scipy.spatial.distance import squareform
from jinja2 import Environment, FileSystemLoader
import re
import warnings
import re

# import tcrdist as td
# from tcrdist.repertoire import TCRrep
# import pandas as pd
# import pwseqdist as pw
# from scipy.spatial.distance import squareform
# from tcrsampler.sampler import TCRsampler
# import palmotif

def dataframe_to_toc_html(df, 
                   cols=['cluster_id','neighbors','K_neighbors'], 
                   template ='toc_template.html',
                   dest='static'):
    """
    Converts a Pandas DataFrame to a toc.html page with links 
    i.e., <a target="logo_frame" href="static/X_cluster_id.html">X</a></li>

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing informaton on each cluster
    cols : list
        List of columns to include in the html version of the dataframe
    template : str
        name of the template in /template/ directory where table is inserted
    dest : str
        name of the directory where the .html will be saved     
    """
    # select the relevant columns
    df = df[cols].copy()

    # add url
    df['link'] = df['cluster_id'].apply(lambda x : f'{x}_cluster_id.html')

    # write html table from DataFrame
    html_table_str = df.to_html(index = False)

    # format the final column to contain a hyperlink
    html_table_str  = _add_format_link(html_str = html_table_str)
    # tbl_body extract the body of the table
    tbl_body = html_table_str.partition("<tbody>\n")[2].partition("</tbody>")[0]
    tbl_body = f"<tbody>\n{tbl_body}</tbody>"

    # tbl_head new head with correct onclick sorts
    tbl_head = _write_th_lines(df)
    
    jinja_env = Environment(loader=FileSystemLoader('templates'))
    table_tmp = jinja_env.get_template(template)
    # insert thead and tbody into template
    html = table_tmp.render(tbl_body = tbl_body, tbl_head = tbl_head)

    # write html to file in destination folder
    name = f'toc.html'
    _write_html(html = html, name = name , dest = dest)

def cluster_dataframe_to_cluster_html(  df, 
                                        cluster_id, 
                                        cols = None, 
                                        template = 'X_cluster_id.html',
                                        dest = 'static'):
    """ ns_df --> static/X_cluster_id.html"""
    
    # select the relevant columns
    if cols is not None:
        df = df[cols].copy()
    
    html_table_str = df.to_html(index = False)
    # tbl_body extract the body of the table
    tbl_body = html_table_str.partition("<tbody>\n")[2].partition("</tbody>")[0]
    tbl_body = f"<tbody>\n{tbl_body}</tbody>"
    # tbl_head new head with correct onclick sorts
    tbl_head = _write_th_lines(df)
    jinja_env = Environment(loader=FileSystemLoader('templates'))
    table_tmp = jinja_env.get_template(template)
    html = table_tmp.render(tbl_body = tbl_body, tbl_head = tbl_head, cluster_id = cluster_id)
    name = f'{cluster_id}_cluster_id.html'
    _write_html(html = html, name = name , dest = dest)

def _write_th_lines(df):
    """
    Write custom header block of html table with correct onclick sort types

    Parameters
    ---------
    df : pd.DataFrame
    
    Returns
    -------
    header : str

    Example
    -------
    >>> _write_th_line(df = pd.DataFrame({"Apple":['mac','ipod'], "Microsoft" : [95,98]}))
    '<th onclick="sortTable(0)">Apple</th>\n<th onclick="sortTableNumeric(1)">Microsoft</th>'
    """
    # get type of first element in each series
    df_dtype = df.apply(lambda x: type(x.iloc[0])).to_list()
    # define a dictionary of sort
    sort_types = {str: 'sortTable', int: 'sortTableNumeric', float:'sortTableNumeric'}
    sort_commands = [sort_types[x] if x in sort_types.keys() else "sortTable" for x in df_dtype]
    th = list()

    for i,name in enumerate(df.columns):
        sort_command = sort_commands[i]
        th.append(f'<th onclick="{sort_command}({i})">{name}</th>')

    header = "\n".join(map(str, th))
    
    return header

def _add_format_link(html_str):
    """ format last column in html table with correct url link type """
    x = html_str.split("\n")
    i = 0 
    edited_html = list()
    for line in x:
        if line.find(".html") != -1:
            gs = re.match(pattern = '(.*<td>)(.*)(</td>)', string = line).groups()
            gs = list(gs)
            link = gs[1]
            #link = line.partition('<td>')[2].partition('</td')[0]
            name = link.split("_cluster_id.html")[0]
            formatted_link = _format_link(link, name = name )
            gs = re.match(pattern = '(.*<td>)(.*)(</td>)', string = line).groups()
            gs = list(gs)
            gs[1] = formatted_link
            line = "".join(gs)
        edited_html.append(line)
    edited_html = "\n".join(edited_html)
    return edited_html 

def _format_link(link, name = None):
    """
    Parameters 
    ----------
    html_str : str
        e.g., to
    Returns
    -------
        s : str
        e.g. 'td><a target="logo_frame" href="LINK">NAME</a></li></td>
    Examples
    --------
    >>> _format_link('static/bkgd_svg.svg')
    '<a target="logo_frame" href="static/bkgd_svg.svg">static/bkgd_svg.svg</a></li>'
    >>> _format_link('static/bkgd_svg.svg', name = 'bkdg')
    '<a target="logo_frame" href="static/bkgd_svg.svg">bkgd</a></li>'
    """
    if name is None:
        name = str(link)
    s = f'<a target="logo_frame" href="{link}">{name}</a></li>'
    return s

# def _write_cluster_table(df, cluster_id,  cols = None, chain = 'beta',dest='static'):
    
#     gene_names = {  'alpha':  ['cdr3_a_aa', 'v_a_gene', 'j_a_gene'],
#                     'beta':   ['cdr3_b_aa', 'v_b_gene', 'j_b_gene'],
#                     'gamma':  ['cdr3_g_aa', 'v_g_gene', 'j_g_gene'],
#                     'delta':  ['cdr3_d_aa', 'v_d_gene', 'j_d_gene']}
    
#     chain_specific_columns = gene_names[chain]
    
#     if cols is None:
#         cols= list()

#     cols =  chain_specific_columns + cols
    
#     # write html table from DataFrame
#     html_table_str = df[cols].to_html(index = False)
#     # extract the body of the table
#     tbody = html_table_str.partition("<tbody>\n")[2].partition("</tbody>")[0]
#     tbody = f"<tbody>\n{tbody}</tbody>"

#     # write a new head with correct onclick sorts
#     thead = _write_th_line(df[cols])
    
#     jinja_env = Environment(loader=FileSystemLoader('templates'))
#     table_tmp = jinja_env.get_template('cluster_table_template.html')
#     # insert thead and tbody into template
#     html = table_tmp.render(tbody = tbody, thead = thead )
#     # write html to file in destination folder
#     name = f'{cluster_id}.df.html'
#     _write_html(html = html, name = name , dest = dest)



def _get_v_j_cdr3( df, cols= None, chain = 'beta'):
    gene_names = {  'alpha':  ['v_a_gene', 'j_a_gene', 'cdr3_a_aa'],
                    'beta':   ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'],
                    'gamma':  ['v_g_gene', 'j_g_gene', 'cdr3_g_aa'],
                    'delta':  ['v_d_gene', 'j_d_gene', 'cdr3_d_aa'] }

    if cols is None:
        cols = gene_names[chain]

    return df[cols]

def _get_cdr3( df, col = None, chain = 'beta'):

    gene_names = {  'alpha': 'cdr3_a_aa',
                    'beta' : 'cdr3_b_aa',
                    'gamma': 'cdr3_g_aa',
                    'delta': 'cdr3_d_aa'}

    if col is None:
        col = gene_names[chain]
    return df[ col ].to_list()

def _get_gene_usage(df, cols = None, chain = 'beta'):

    gene_names = {  'alpha': ['v_a_gene', 'j_a_gene'],
                    'beta' : ['v_b_gene', 'j_b_gene'],
                    'gamma': ['v_g_gene', 'j_g_gene'],
                    'delta': ['v_d_gene', 'j_d_gene']}
    if cols is None:
        cols = gene_names[chain]
    
    gene_usage = df.groupby( cols ).size()

    sample_usage = gene_usage.reset_index().to_dict('split')['data']

    return sample_usage

def _strip_allele(gene_string):
    """
    Example
    -------
    >>> _strip_allele('TRBV12-2*01'
    >>> 'TRBV12-2'

    >>> df[['v','j']] = df[['v','j']].apply(lambda x: x.apply(_strip_allele))
    """
    try:
        regex_groups = re.match(string = gene_string , pattern ='(.*)\*')
        return regex_groups[1]
    except AttributeError:
        return gene_string

def _pick_best(gene_string):
    """
    Example
    -------
    >>> _pick_best('TRAV12-2*01,TRAV12D-1*01')
    TRAV12-2*01
    >>>_pick_best('TRAV12-2*01'
    'TRAV12-2*01'
    """
    return gene_string.split(",")[0]

def _add_allele(gene_string, allele = "*01"):
    """
    Example
    -------
    >>> _add_allele('TRBV12-2')
    'TRBV12-2*01'
    """
    return f"{gene_string}{allele}"


def get_centroid_seq(seqs, metric = pw.metrics.nw_hamming_metric):
    """
    Given a list of sequences, returns the sequence with the minimum 
    sum of distances to all other seqs in the list.

    Parameters
    ----------
    seqs : list
        list of strings (amino acid rep)
    metric : func
        defaults to pw.metrics.nw_hamming_metric

    Returns
    -------
    centroid_seq : str

    Example 
    -------
    >>> seqs = ['CASSEILAALGTQYF', 'CASSWTSRETQYF', 'CASSLAQETQYF', 'CASSLAPGDVSQYF', 'CASSWDQETQYF', 'CASSLWWDSGANVLTF', 'CASSLARTLSSGANVLTF', 'CASIPGTLFTFSGANVLTF', 'CASSFASSGANVLTF', 'CASSYRLLSGANVLTF']	
    >>> get_centroid_seq(seqs)
    'CASSFASSGANVLTF'

    Notes 
    -----
    In case of multiple occurrences of the minimum values, the indices 
    corresponding to the first occurrence are returned.
    """
    #import pwseqdist as pw
    #from scipy.spatial.distance import squareform
    if len(seqs) < 3:
        return seqs[0]
    dvec = pw.apply_pairwise_sq(seqs = seqs, 
                metric = metric , 
                ncpus  = 1 )
    dmat = squareform(dvec).astype(int)
    index_pos_of_min = dmat.sum(axis = 0).argmin()
    centroid_seq  = seqs[index_pos_of_min]
    return centroid_seq

def _get_neighbors(df, ids):
    """ 
    Returns a subset of multiple rows from a 
    DataFame < df > as a copy < df_subset >. 
    The rows are defined by list of indices < ids >.

    Parameters
    ----------
    df : DataFrame
    
    idx : list

    Returns
    -------
    df_subset : DataFrame
        subset of input < df >

    Notes
    -----
    WARNING: This only works if 'neighbor' column produces a list of indices 
    matching clone_df indices 
    """
    df_subset = df.iloc[ids,:].copy()
    return df_subset

def _get_neighbors_from_row(row, df, col = 'neighbors'):
    """ 
    Returns a subset of multiple rows from a DataFame < df > as a copy < df_subset >.

    Uses an the attriute, specified by < col > , in a Pandas Series < row >,
    containing a lists of indices used to define the 
    subset of < df >.

    Parameters
    ----------
    row : row of pd.DataFrame
    df : pd.DataFrame
    
    Returns
    -------
    df_subset : pd.DataFrame
        subset of input < df >
    
    WARNING: This only works if 'neighbor' column produces a list of indices 
    matching clone_df indices 
    """
    idx = row[col]
    if isinstance(idx, str):
        idx = list(pd.eval(idx)) 

    df_subset = _get_neighbors(df, idx)
    return df_subset

def _get_cluster_attributes(df, 
                            row= None, 
                            ids = None, 
                            col = 'neighbors', 
                            chain = 'beta',
                            centroid_metric = pw.metrics.nw_hamming_metric):
    """
    Prepare a named tuple of all needed elements for running palmotif on a cluster of clones. 

    Parameters 
    ----------
    df : pd.DataFrame
        TCRrep.clone_df of a Pandas DataFrame containing clone_df type information
    row : None or pd.Series
        Pandas Sereies containg column 
    ids : None or list 
        List of indices to select from < df >. e.g., [1, 4, 10, 100]
    col : string or None
        If row is not None, this specifies column containing list of cluster neighbor or 
        members. If row is None, this argument is not used.
    chain : str
        'beta', 'alpha', 'gamma' or 'delta'
    centroid_metric
    
    Returns
    -------
    cl : namedtuple
        cluster_id : int or None
        neighbors : list
        K_neighbors : int
        ns_v_j_cdr3 : pd.DataFrame
        ns_df :pd.DataFrame
        ns_cdr3_seqs : list e.g., ['CASSALTGPAIAKNIQYF', 'CASSKLTGEAVAKNIQYF', 'CASSPLTGEALAKNIQYF'] 
        gene_usage : list of lists e.g., [['TRBV27*01', 'TRBJ2-4*01', 9]]
        centroid : str 
    
    Notes
    -----
    Any clustering method can be used and indices provided directly to ids, but these must match
    indices in clone_df.
    """
    if row is not None:
        ns_df     = _get_neighbors_from_row(row = row, df = df, col = col)
        cluster_id = row['cluster_id']
    else:
        ns_df    = _get_neighbors(ids = ids, df = df)
        cluster_id = None
    #<ns_v_j_cdr3 > neighbors v j cdr3
    ns_v_j_cdr3  = _get_v_j_cdr3(df = ns_df, chain = chain)
    # <ns_seqs> neighbors cdr3 seqs
    ns_cdr3_seqs = _get_cdr3(df = ns_df, chain = chain) 
    # <ns_seqs> Gene Usage in Format for sampling (e.g., [['TRBV27*01', 'TRBJ2-4*01', 9]])
    gene_usage_raw   = _get_gene_usage(df = ns_v_j_cdr3, chain = chain)
    # change all alleles to *01
    gene_usage       = [[_strip_allele(x[0]), _strip_allele(x[1]), x[2]] for x in gene_usage_raw]
    gene_usage       = [[_add_allele(x[0],  allele = "*01"), _add_allele(x[1], allele = "*01"), x[2]] for x in gene_usage]
    
    # <centroid> CDR3 with smallest total distance to other CDR3s, used for motif alignment
    centroid     = get_centroid_seq(seqs=ns_cdr3_seqs, metric = centroid_metric)
    
    cluster_nt = namedtuple('cluster', ['cluster_id','neighbors','K_neighbors','ns_df','ns_v_j_cdr3','ns_cdr3_seqs','gene_usage','centroid'])
    
    neighbors = ns_df.index.to_list()
    K_neighbors = len(neighbors)

    cl = cluster_nt(cluster_id = cluster_id,
                 neighbors = neighbors,
                 K_neighbors =  K_neighbors,
                 ns_df = ns_df , 
                 ns_v_j_cdr3 = ns_v_j_cdr3, 
                 ns_cdr3_seqs = ns_cdr3_seqs, 
                 gene_usage = gene_usage, 
                 centroid = centroid)
    return cl

def _sample(gene_usage, ts, depth = 100):
    """
    gene_usage : list
    list of lists 
    ts : tcrsampler 
    """
    with warnings.catch_warnings(record=True) as w:
        x = ts.sample(gene_usage, depth = depth, flatten = True)

    try: 
        w = str(w[-1].message)
    except:
        w = None
    
    sample = [i for i in x if i is not None] 
    return sample, w  

def _modify_svg(svg, w ='1000px'):
    svg = svg.replace('height="100%"', f'width="{w}"')
    svg = svg.replace('width="100%"', "")
    return svg  

class MotifReport():
    def __init__(self, cluster_attributes):
        """ Class for holding and writing cluster motif/svg objects """
        self.cluster_attributes = cluster_attributes
        self.warning     = None
        self.raw_motif   = None
        self.raw_stat    = None
        
        self.pal_motif   = None
        self.pal_stat    = None
        
        self.bkgd_motif  = None
        self.bkgd_stat   = None

        self.raw_svg  = None
        self.pal_svg  = None
        self.bkgd_svg = None

    def _write_svg(self, svg, name, dest = 'static'):
        if not os.path.isdir(dest):
            os.mkdir(dest)

        fn = os.path.join(dest, name)
        with open(fn, 'w') as fh:
            fh.write(f"<body>{svg}</body>")
    
    def _write_svgs(self, svgs, name, dest = 'static'):
        import os
        fn = os.path.join(dest, name)
        with open(fn, 'w') as fh:
            fh.write(f"<body>\n")
            for svg in svgs:
                fh.write(f"\t{svg}\n")
            fh.write(f"\n</body>")   


    def _write_palmotifs(self, dest = 'static'):
        cluster_id = self.cluster_attributes.cluster_id
        
        raw_name = f"{cluster_id}.raw_svg.svg"
        self._write_svg(self.raw_svg, name = raw_name, dest=dest)
        
        pal_name = f"{cluster_id}.pal_svg.svg"
        self._write_svg(self.pal_svg, name = pal_name, dest=dest)
        
        bkgd_name = f"{cluster_id}.bkgd_svg.svg"
        self._write_svg(self.bkgd_svg, name = bkgd_name, dest=dest)
    
    def _write_motif_html(self, dest = 'static'):

        jinja_env = Environment(loader=FileSystemLoader('templates'))
        # the row template is for the motifs
        row_tmp = jinja_env.get_template('row_template.html')
        # the page template is for the 
        page_tmp = jinja_env.get_template('2col_template.html')
    
        row_html = row_tmp.render(x= f"rank{rank}.clone{i}")
        
        #page_html = page_tmp.render(x = row_html)
        #_write_html(html = page_html, name = f"rank{rank}.clone{i}.view.html", dest =".")
        

def _get_palmotifs(cluster_attributes, sampler, write = True, dest = 'static3'):
    """
    Use cluster attributes to generate palmotif matrices and svg
    """
    
    cl = cluster_attributes
    
    refs, w = _sample(gene_usage = cl.gene_usage, ts = sampler)
    
    bkgd = list(refs) + [cl.centroid] 

    m = MotifReport(cluster_attributes = cl)
    m.warning = w
    m.raw_motif  = palmotif.compute_pal_motif(seqs = cl.ns_cdr3_seqs, centroid = cl.centroid)
    m.pal_motif  = palmotif.compute_pal_motif(seqs = cl.ns_cdr3_seqs, refs= refs, centroid = cl.centroid)
    m.bkgd_motif = palmotif.compute_pal_motif(seqs = bkgd , centroid = cl.centroid)

    m.raw_svg    = palmotif.svg_logo(m.raw_motif[0],  return_html=False, return_str = True, svg_height='100px', svg_width='800px')
    m.pal_svg    = palmotif.svg_logo(m.pal_motif[0],  return_html=False, return_str = True, svg_height='100px', svg_width='800px')
    m.bkgd_svg   = palmotif.svg_logo(m.bkgd_motif[0], return_html=False, return_str = True, svg_height='100px', svg_width='800px')

    m.raw_svg    = _modify_svg(m.raw_svg,  w ='250px height=50px')
    m.pal_svg    = _modify_svg(m.pal_svg,  w ='250px height=50px')
    m.bkgd_svg   = _modify_svg(m.bkgd_svg, w ='200px height=50px')

    m.raw_stat   = m.raw_motif[1]
    m.pal_stat   = m.pal_motif[1]
    m.bkgd_stat  = m.bkgd_motif[1]
    
    if write:
        #m._write_palmotifs(dest = dest)
        name = f'{m.cluster_attributes.cluster_id}_cluster_id.svg'
        m._write_svgs(svgs = [m.pal_svg, m.raw_svg, m.bkgd_svg], name=name, dest = 'static')

    return m


def _write_html(html, name, dest = 'static3'):
    if not os.path.isdir(dest):
        os.mkdir(dest)
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write(f"{html}")

def _write_svgs(svgs, name, dest = 'static'):
	import os
	fn = os.path.join(dest, name)
	with open(fn, 'w') as fh:
		fh.write(f"<body>\n")
		for svg in svgs:
			fh.write(f"\t{svg}\n")
		fh.write(f"\n</body>")   


def _write_svgs_w_labs(svgs, labs, name, dest = 'static'):
    import os
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write(f"<body>\n")
        for svg,lab in zip(svgs, labs):
            fh.write(f"\t{svg}\n")
            fh.write(f"<div><p>{lab}</p></div>\n")
        fh.write(f"\n</body>")   









# import tcrdist as td
# from tcrdist.repertoire import TCRrep
# import pandas as pd
# import pwseqdist as pw
# from scipy.spatial.distance import squareform
# from tcrsampler.sampler import TCRsampler
# import palmotif
# from collections import namedtuple
# import os



# def aggregate_cluster_graphics( tr,
#                                 min_K_neighors = 5, 
#                                 chain = 'beta',
#                                 centroid_metric = pw.metrics.nw_hamming_metric,
#                                 ):

#     clusters = list()
#     for i,row in tr.cluster_df.iterrows():
#         if row.K_neighbors > min_K_neighors:
#             cl = cluster._get_cluster_attributes(df = tr.clone_df, 
#                                                 row = row,
#                                                 chain = chain,
#                                                 metric = centroid_metric)
            
#             cm = cluster._tcrdist2_motif(cluster_attributes = cl, sampler = t)
#             clusters.append(cm)
    
#     return clusters

# def _write_svg(svg, name, dest = 'static'):
#     fn = os.path.join(dest, name)
#     with open(fn, 'w') as fh:
#         fh.write(f"<body>{svg}</body>")

# def _write_html(html, name, dest = 'static'):
#     fn = os.path.join(dest, name)
#     with open(fn, 'w') as fh:
#         fh.write(f"{html}")

# def write_cluster_graphics(clusters, dest = "static3"):
#     if not os.path.isdir(dest):
#         os.mkdir(dest)
#     for cl in clusters:
#         for n in ['raw_svg','pal_svg', 'bkgd_svg']:
#             id = cl.cluster_attributes.cluster_id
#             svg_name = f"{id}.{name}

#             _write_svg(cl.__getattribute__(n), name = n, dest = dest)



    



# class MotifReport():
#     def __init__(self, cluster_attributes):
#         """ Class for holding cluster motif objects """
#         self.cluster_attributes = cluster_attributes
        
#         self.raw_motif   = None
#         self.raw_likl    = None
#         self.pal_motif   = None
#         self.pal_likl    = None
#         self.bkgd_motif  = None
#         self.bkgd_likl   = None

#         self.raw_svg  = None
#         self.pal_svg  = None
#         self.bkgd_svg = None


# def _get_cluster_attributes(df, 
#                             row= None, 
#                             ids = None, 
#                             col = 'neighbors', 
#                             chain = 'beta',
#                             metric = pw.metrics.nw_hamming_metric):
#     """
#     Prepare a named tuple of all needed 
#     elements for running palmotif on a cluster of clones. 

#     Parameters 
#     ----------
#     df : pd.DataFrame
#         TCRrep.clone_df of a Pandas DataFrame containing clone_df type information
#     row : None or pd.Series
#         Pandas Sereies containg column 
#     ids : None or list 
#         List of indices to select from < df >. e.g., [1, 4, 10, 100]
#     col : string or None
#         If row is not None, this specifies column containing list of cluster neighbor or 
#         members. If row is None, this argument is not used.
#     chain = 'beta'
    
#     Return
#     ------
#     cl : namedtuple
#         ns_v_j_cdr3 : pd.DataFrame
#         ns_df :pd.DataFrame
#         ns_cdr3_seqs : list e.g., ['CASSALTGPAIAKNIQYF', 'CASSKLTGEAVAKNIQYF', 'CASSPLTGEALAKNIQYF'] 
#         gene_usage : list of lists e.g., [['TRBV27*01', 'TRBJ2-4*01', 9]]
#         centroid : str 
    
#     Notes
#     -----
#     Any clustering method can be used and indices provided directly to ids, but these must match
#     indices in clone_df or metadata df.
#     """
#     if row is not None:
#         ns_df     = _get_neighbors_from_row(row = row, df = df, col = col)
#         cluster_id = row['cluster_id']
#     else:
#         ns_df    = _get_neighbors(ids = ids, df = df)
#         cluster_id = None
#     #<ns_v_j_cdr3 > neighbors v j cdr3
#     ns_v_j_cdr3  = _get_v_j_cdr3(df = ns_df, chain = chain)
#     # <ns_seqs> neighbors cdr3 seqs
#     ns_cdr3_seqs = _get_cdr3(df = ns_df, chain = chain) 
#     # <ns_seqs> Gene Usage in Format for sampling (e.g., [['TRBV27*01', 'TRBJ2-4*01', 9]])
#     gene_usage   = _get_gene_usage(df = ns_v_j_cdr3, chain = chain)
#     # <centroid> CDR3 with smallest total distance to other CDR3s, used for motif alignment
#     centroid     = get_centroid_seq(seqs=ns_cdr3_seqs, metric = metric)
    
#     cluster = namedtuple('cluster', ['neighbors','K_neighbors','ns_df','ns_v_j_cdr3','ns_cdr3_seqs','gene_usage','centroid'])
    
#     neighbors = ns_df.index.to_list()
#     K_neighbors = len(neighbors)

#     cl = cluster(cluster_id = cluster_id,
#                  neighbors = neighbors,
#                  K_neighbors =  K_neighbors,
#                  ns_df = ns_df , 
#                  ns_v_j_cdr3 = ns_v_j_cdr3, 
#                  ns_cdr3_seqs = ns_cdr3_seqs, 
#                  gene_usage = gene_usage, 
#                  centroid = centroid)
#     return cl



# def _get_v_j_cdr3( df, cols= None, chain = 'beta'):
#     gene_names = {  'alpha':  ['v_a_gene', 'j_a_gene', 'cdr3_a_aa'],
#                     'beta':   ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'],
#                     'gamma':  ['v_g_gene', 'j_g_gene', 'cdr3_g_aa'],
#                     'delta':  ['v_d_gene', 'j_d_gene', 'cdr3_d_aa'] }

#     if cols is None:
#         cols = gene_names[chain]

#     return df[cols]

# def _get_cdr3( df, col = None, chain = 'beta'):

#     gene_names = {  'alpha': 'cdr3_a_aa',
#                     'beta' : 'cdr3_b_aa',
#                     'gamma': 'cdr3_g_aa',
#                     'delta': 'cdr3_d_aa'}

#     if col is None:
#         col = gene_names[chain]
#     return df[ col ].to_list()

# def _get_gene_usage(df, cols = None, chain = 'beta'):

#     gene_names = {  'alpha': ['v_a_gene', 'j_a_gene'],
#                     'beta' : ['v_b_gene', 'j_b_gene'],
#                     'gamma': ['v_g_gene', 'j_g_gene'],
#                     'delta': ['v_d_gene', 'j_d_gene']}
#     if cols is None:
#         cols = gene_names[chain]
    
#     gene_usage = df.groupby( cols ).size()

#     sample_usage = gene_usage.reset_index().to_dict('split')['data']

#     return sample_usage

# def get_centroid_seq(seqs, metric = pw.metrics.nw_hamming_metric):
#     """
#     Given a list of sequences, return the sequences with the minimum sum of distances to all other seqs in the list.

#     Parameters
#     ----------
#     seqs : list
#         list of strings (amino acid rep)
#     metric : func
#         defaults to pw.metrics.nw_hamming_metric

#     Returns
#     -------
#     centroid_seq : str

#     Example 
#     -------
#     >>> seqs = ['CASSEILAALGTQYF', 'CASSWTSRETQYF', 'CASSLAQETQYF', 'CASSLAPGDVSQYF', 'CASSWDQETQYF', 'CASSLWWDSGANVLTF', 'CASSLARTLSSGANVLTF', 'CASIPGTLFTFSGANVLTF', 'CASSFASSGANVLTF', 'CASSYRLLSGANVLTF']	
#     >>> get_centroid_seq(seqs)
#     'CASSFASSGANVLTF'

#     Notes 
#     -----
#     In case of multiple occurrences of the minimum values, the indices corresponding to the first occurrence are returned.
#     """
#     #import pwseqdist as pw
#     #from scipy.spatial.distance import squareform
#     dvec = pw.apply_pairwise_sq(seqs = seqs, 
#                 metric = metric , 
#                 ncpus  = 1 )
#     dmat = squareform(dvec).astype(int)
#     index_pos_of_min = dmat.sum(axis = 0).argmin()
#     centroid_seq  = seqs[index_pos_of_min]
#     return centroid_seq

    





     

