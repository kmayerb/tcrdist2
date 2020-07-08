""" 
cluster.py

Functions for clustering TCR clones and MotifReport STORAGE class for holding cluster specific 
attributes for later recall.

e.g., motifs is dictionary containing MotifReport instances
{51: <tcrdist.cluster.MotifReport at 0x1356ac610>,
 68: <tcrdist.cluster.MotifReport at 0x137f05750>,
 50: <tcrdist.cluster.MotifReport at 0x134c370d0>}

motifs[51].pal_motif # accesses palmotif tuple)
motifs[51].pal_stat  # access kl and loglik stats
motifs[51].cluster_attributes # access a namedtuple set of attributes used elsewhere when making a report

"""

from collections import namedtuple
from jinja2 import Environment, FileSystemLoader, PackageLoader
import os
import palmotif
import pandas as pd
import pwseqdist as pw
from scipy.spatial.distance import squareform
import re
import warnings

class MotifReport():
    def __init__(self, cluster_attributes):
        """ 
        Simple class for holding a number of cluster motif/svg objects
        
        Attributes
        ----------
        cluster_attributes : namedtuple 
            namedtuple generated by _get_cluster_attributes() 
        raw_motif : palmotif tuple
            palmotif WITHOUT background specified generated by 
            palmotif.compute_pal_motif in _get_palmotifs() 
        pal_motif : palmotif tuple
            palmotif WITH background generated by 
            palmotif.compute_pal_motif in  _get_palmotifs() 
        bkgd_motif : pamotif tuple
            palmotif using the centroid and background seqs only 
            generated by palmotif.compute_pal_motif in  _get_palmotifs() 
        pal_stat : pd.DataFrame 
            kl and loglik
        warning : str
            holds warning message potentially captured from tcrsampler

        Notes 
        -----
        cluster_attributes has the following attributes: 
            cluster_id : int or None
            neighbors : list
            K_neighbors : int
            ns_v_j_cdr3 : pd.DataFrame
            ns_df :pd.DataFrame
            ns_cdr3_seqs : list e.g., ['CASSALTGPAIAKNIQYF', 'CASSKLTGEAVAKNIQYF', 'CASSPLTGEALAKNIQYF'] 
            gene_usage_raw : list of lists e.g., [['TRBV27*02', 'TRBJ2-4*01', 9]]
            gene_usage : list of lists e.g., [['TRBV27*01', 'TRBJ2-4*01', 9]]
            centroid : str 

        """
        # namedtuple of cluster attributes
        self.cluster_attributes = cluster_attributes
    
        # Motifs
        self.raw_motif   = None
        self.raw_stat    = None
        self.pal_motif   = None
        self.pal_stat    = None
        self.bkgd_motif  = None
        self.bkgd_stat   = None
        # SVG attributes
        self.raw_svg  = None
        self.pal_svg  = None
        self.bkgd_svg = None
        #Warning
        self.warning     = None

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
        Has the following attributes: 
            cluster_id : int or None
            neighbors : list
            K_neighbors : int
            ns_v_j_cdr3 : pd.DataFrame
            ns_df :pd.DataFrame
            ns_cdr3_seqs : list e.g., ['CASSALTGPAIAKNIQYF', 'CASSKLTGEAVAKNIQYF', 'CASSPLTGEALAKNIQYF'] 
            gene_usage_raw : list of lists e.g., [['TRBV27*02', 'TRBJ2-4*01', 9]]
            gene_usage : list of lists e.g., [['TRBV27*01', 'TRBJ2-4*01', 9]]
            centroid : str 
        
    Notes
    -----
    Any clustering method can be used and indices provided directly to ids, but these must match
    indices in clone_df.

    WARNING: THE DEFAULT BEHAIVIOR IS TO SET ALL ALLELES TO *01. in gene_usage. Original alleles are 
    in gene_usage_raw attribute. Because most high throughput datasets lack specific allele information
    tcrsampler typically repressents all clones as coming from the *01 allele.

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
                 neighbors     = neighbors,
                 K_neighbors   =  K_neighbors,
                 ns_df         = ns_df , 
                 ns_v_j_cdr3   = ns_v_j_cdr3, 
                 ns_cdr3_seqs  = ns_cdr3_seqs, 
                 gene_usage    = gene_usage, 
                 centroid      = centroid)
    return cl

def _get_palmotifs(cluster_attributes, sampler, depth = 1000, write = True, dest = 'static'):
    """
    Use cluster attributes to generate palmotif matrices and svg
    """
    
    cl = cluster_attributes
    
    refs, w = _sample(gene_usage = cl.gene_usage, ts = sampler, depth = depth)
    
    bkgd = list(refs) + [cl.centroid] 

    m = MotifReport(cluster_attributes = cl)
    m.warning = w
    m.raw_motif  = palmotif.compute_pal_motif(seqs = cl.ns_cdr3_seqs, centroid = cl.centroid)
    m.pal_motif  = palmotif.compute_pal_motif(seqs = cl.ns_cdr3_seqs, refs= refs, centroid = cl.centroid)
    m.bkgd_motif = palmotif.compute_pal_motif(seqs = bkgd , centroid = cl.centroid)

    m.raw_svg    = palmotif.svg_logo(m.raw_motif[0],  return_html=False, return_str = True, svg_height='100px', svg_width='800px')
    m.pal_svg    = palmotif.svg_logo(m.pal_motif[0],  return_html=False, return_str = True, svg_height='100px', svg_width='800px')
    m.bkgd_svg   = palmotif.svg_logo(m.bkgd_motif[0], return_html=False, return_str = True, svg_height='100px', svg_width='800px')

    m.raw_svg    = _modify_svg(m.raw_svg,  w ='150px height=40px')
    m.pal_svg    = _modify_svg(m.pal_svg,  w ='150px height=40px')
    m.bkgd_svg   = _modify_svg(m.bkgd_svg, w ='100px height=40px')

    m.raw_stat   = m.raw_motif[1]
    m.pal_stat   = m.pal_motif[1]
    m.bkgd_stat  = m.bkgd_motif[1]
    
    if write:
        #m._write_palmotifs(dest = dest)
        name = f'{m.cluster_attributes.cluster_id}_cluster_id.svg'
        m._write_svgs(svgs = [m.pal_svg, m.raw_svg, m.bkgd_svg], name=name, dest = 'static')

    return m





def dataframe_to_toc_html(df, 
                   cols=['cluster_id','neighbors','K_neighbors'], 
                   name = 'toc.html',
                   template ='toc_template.html',
                   dest='static',
                   write = True,
                   jinja_local = False):
    """
    Converts Pandas DataFrame, with specified columns,
    to a toc.html page. 

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing informaton on each cluster
    cols : list
        List of columns to include in the html version of the dataframe
    name : str
        Name of the output file (by default /static/toc.html)
    template : str
        name of the template in /template/ directory where table is inserted
    dest : str
        name of the directory where the .html will be saved
    write : bool 
        If True writes file < name > to < dest > folder 
    jinja_local : bool
        Set to True if you want ot use templates in working dir. Set to False to 
        use template in the tcrdist2 package

    Returns 
    -------
    html : str
        html page   

    Examples
    --------
    >>> from tcrdist import cluster
    >>> html = cluster.dataframe_to_toc_html(df = TCRrep.cluster_df)
    
    Notes
    -----
    The page contains links 
    i.e., <a target="logo_frame" href="static/x_cluster_id.html">X</a></li> 

    You may have to set jinja_local to True for CI testing.
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
    
    if jinja_local:
        # if you want to use a local directory called templates
        jinja_env = Environment(loader=FileSystemLoader('templates'))
    else:
        # if you want to use templates within the tcrdist2 package
        jinja_env = Environment(loader=PackageLoader('tcrdist', 'templates'))
    
    table_tmp = jinja_env.get_template(template)
    # insert thead and tbody into template
    html = table_tmp.render(tbl_body = tbl_body, tbl_head = tbl_head)

    if write:
        # write html to file in destination folder
        _write_html(html = html, name = name , dest = dest)
    
    return html

def cluster_dataframe_to_cluster_html(  df, 
                                        cluster_id, 
                                        cols = None, 
                                        template = 'X_cluster_id.html',
                                        dest = 'static',
                                        write = True,
                                        jinja_local = False):
    """ 
    Converts Pandas DataFrame with cluster-specific information 
    into an html page. Output name is set to f'{cluster_id}_cluster_id.html'
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with cluster specific information 
    cluster_id : int
        Output name is set to f'{cluster_id}_cluster_id.html'
    cols : None or list
        List of colums to present in the html table. If None, all columns are used.
    template : str
        Name of template in the jinja templates folder.
    write : bool
        If True writes file < name > to < dest > folder
    jinja_local : bool
        Set to True if you want ot use templates in working dir. Set to False to 
        use template in the tcrdist2 package
    
    Returns
    -------
    html : str
        html page from template incoporating table
    
    Notes
    -----
    https://jinja.palletsprojects.com/en/2.11.x/api/#loaders

    You may have to set jinja_local to True for CI testing.
    """

    # select the relevant columns
    if cols is not None:
        df = df[cols].copy()
    
    html_table_str = df.to_html(index = False)
    # tbl_body extract the body of the table
    tbl_body = html_table_str.partition("<tbody>\n")[2].partition("</tbody>")[0]
    tbl_body = f"<tbody>\n{tbl_body}</tbody>"
    # tbl_head new head with correct onclick sorts
    tbl_head = _write_th_lines(df)

    if jinja_local:
        # if you want to use a local directory called templates
        jinja_env = Environment(loader=FileSystemLoader('templates'))
    else:
        # if you want to use templates within the tcrdist2 package
        jinja_env = Environment(loader=PackageLoader('tcrdist', 'templates'))

    table_tmp = jinja_env.get_template(template)
    html = table_tmp.render(tbl_body = tbl_body, tbl_head = tbl_head, cluster_id = cluster_id)
    name = f'{cluster_id}_cluster_id.html'
    if write: 
        _write_html(html = html, name = name , dest = dest)
    return html

def _write_th_lines(df):
    """
    Write custom <th>header</th> block of html table with correct onclick sort types. 
    Numeric columns need sortTableNumeric(i) and alpha numeric need sortTable(i). 

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
    """ 
    Adds formatted links the last column in an html table. 
    Each link is added by function _format_link (see below)
    
    Parameter
    ---------
    html_str : str

    Returns
    -------
    edited_html : str
    """
    x = html_str.split("\n")
    i = 0 
    edited_html = list()
    for line in x:
        if line.find(".html") != -1:
            gs = re.match(pattern = r'(.*<td>)(.*)(</td>)', string = line).groups()
            gs = list(gs)
            link = gs[1]
            #link = line.partition('<td>')[2].partition('</td')[0]
            name = link.split("_cluster_id.html")[0]
            formatted_link = _format_link(link, name = name )
            gs = re.match(pattern = r'(.*<td>)(.*)(</td>)', string = line).groups()
            gs = list(gs)
            gs[1] = formatted_link
            line = "".join(gs)
        edited_html.append(line)
    edited_html = "\n".join(edited_html)
    return edited_html 

def _format_link(link, name = None):
    """
    Formats links destination to a hyperlink (see example)

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

def _get_v_j_cdr3( df, cols= None, chain = 'beta'):
    """ 
    Parameters
    ----------
    df : pd.DataFame
        Data frame containing cdr3 seqeunce and v and j gene names
    cols : list or None
        list of columns names, not needed if chain is specified
    chain : str
        Must be one of the following : ['beta', 'alpha', 'gamma', 'delta']

    Returns 
    -------
    df_cols : pd.DataFrame
        Dataframe with specific columns.

    Notes
    ------
    handles chain specific names
    """
    gene_names = {  'alpha':  ['v_a_gene', 'j_a_gene', 'cdr3_a_aa'],
                    'beta':   ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'],
                    'gamma':  ['v_g_gene', 'j_g_gene', 'cdr3_g_aa'],
                    'delta':  ['v_d_gene', 'j_d_gene', 'cdr3_d_aa'] }

    if cols is None:
        cols = gene_names[chain]

    df_cols = df[cols].copy()

    return df_cols

def _get_cdr3( df, col = None, chain = 'beta'):
    """
    Parameters
    ----------
    df : pd.DataFame
        Data frame containing cdr3 seqeunce and v and j gene names
    col: str
        column name containing cdr3 seqs; however, not needed if chain is specified
    chain : str
        Must be one of the following : ['beta', 'alpha', 'gamma', 'delta']

    Returns 
    -------
    cdr3_list : list 
        list of cdr3 strings
    """
    gene_names = {  'alpha': 'cdr3_a_aa',
                    'beta' : 'cdr3_b_aa',
                    'gamma': 'cdr3_g_aa',
                    'delta': 'cdr3_d_aa'}

    if col is None:
        col = gene_names[chain]
    
    cdr3_list = df[ col ].to_list()

    return cdr3_list


def _get_gene_usage(df, cols = None, chain = 'beta'):
    """
    Returns a list of list specifying frequency of V-J gene pairs

    Parameters
    ----------
    df : pd.DataFame
        Data frame containing cdr3 seqeunce and v and j gene names
    cols : list or None
        list of columns names, not needed if chain is specified
    chain : str
        Must be one of the following : ['beta', 'alpha', 'gamma', 'delta']

    Returns 
    -------
    sample usage : list 
        list of lists with pairs of v and j genes and number of times that pairing is observed 
        e.g., [['TRBV13-2*01', 'TRBJ2-3*01', 1]]
    """
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
    Removes allele (see example)

    Example
    -------
    >>> _strip_allele('TRBV12-2*01')
    >>> 'TRBV12-2'

    To do this to many columns of a dataframe:
    >>> df[['v','j']] = df[['v','j']].apply(lambda x: x.apply(_strip_allele))
    """
    try:
        regex_groups = re.match(string = gene_string , pattern =r'(.*)\*')
        return regex_groups[1]
    except AttributeError:
        return gene_string

def _pick_best(gene_string):
    """
    Pick the best (in this case first) gene usage in a comma separated list. 

    Example
    -------
    >>> _pick_best('TRAV12-2*01,TRAV12D-1*01')
    TRAV12-2*01
    >>>_pick_best('TRAV12-2*01')
    'TRAV12-2*01'
    """
    return gene_string.split(",")[0]

def _add_allele(gene_string, allele = "*01"):
    """
    Add an allele number to a gene name (see example)
    
    Example
    -------
    >>> _add_allele('TRBV12-2')
    'TRBV12-2*01'
    """
    return f"{gene_string}{allele}"


def get_centroid_seq(seqs, metric = pw.metrics.nw_metric):
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
    dmat = pw.apply_pairwise_sq(seqs = seqs, 
                metric = metric , 
                ncpus  = 1 )
    dmat = dmat.astype(int)
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

    Parameters
    ----------
    row : row of pd.DataFrame
        
    df : pd.DataFrame
    
    Returns
    -------
    df_subset : pd.DataFrame
        subset of input < df >
    
    Notes
    -----
    Uses an the attriute, specified by < col > , in a Pandas Series < row >,
    containing a lists of indices used to define the 
    subset of < df >.

    WARNING: This only works if 'neighbor' column produces a list of indices 
    matching clone_df indices 

    """
    idx = row[col]
    if isinstance(idx, str):
        idx = list(pd.eval(idx)) 

    df_subset = _get_neighbors(df, idx)
    return df_subset


def _sample(gene_usage, ts, depth = 1000):
    """
    Parameters
    ----------
    gene_usage : list
        list of lists usually accessed from MotifReport.cluster_attributes.gene_usage 
    ts : tcrsampler 
        tcrsampler instance that has been initialized and is ready for sampling
    depth : int 
        For every real sequence in the cluster how many sequences should be sampled to generate a background
    
    Returns 
    -------
    sample : list of cdr3 seqs from background set
    
    w : str
        warning message string from tcrsampler (could be alerting that v,j pair is missing), 
        which we want to be able to propogate to the MotifReport instance

    """
    with warnings.catch_warnings(record=True) as w:
        x = ts.sample(gene_usage, depth = depth, flatten = True)
    try: 
        w = str(w[-1].message)
    except:
        w = None
    
    sample = [i for i in x if i is not None] 
    return sample, w  
        

def _write_html(html, name, dest = 'static'):
    """ 
    writes a html string to a file dest/name 
    
    Parameters 
    ----------
    html : str
        string of html text
    name : str
        filename to output html string
    dest : str
        destination folder
    """

    if not os.path.isdir(dest):
        os.mkdir(dest)
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write(f"{html}")

def _write_svgs(svgs, name, dest = 'static'):
    """ 
    writes a svg string from multiple svgs to a file dest/name 
    
    Parameters 
    ----------
    svgs : list of strs
        list of svgs strings
    name : str
        filename to output svg composite string
    dest : str
        destination folder
    """
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write(f"<body>\n")
        for svg in svgs:
            fh.write(f"\t{svg}\n")
        fh.write(f"\n</body>")   

def _write_svgs_w_labs(svgs, labs, name, dest = 'static', hl = 'h1'):
    """ 
    writes a svg string from multiple svgs and multiple labels to a file dest/name 
    
    Parameters 
    ----------
    svgs : list of strs
        list of svgs strings
    name : str
        filename to output svg composite string
    dest : str
        destination folder
    hl : str
        header level e.g., 'h1' , 'h2', 'h3'
    
    """
    fn = os.path.join(dest, name)
    with open(fn, 'w') as fh:
        fh.write('<!doctype html>\n<html lang="en">')
        fh.write(f'<{hl}>Cluster Tree</{hl}></div>\n')
        fh.write(f"<body>\n")

        for svg,lab in zip(svgs, labs):
            fh.write(f'<div><p class="caption">{lab}</p></div>\n')
            fh.write(f"\t{svg}\n")
            
        fh.write(f"\n</body>\n") 
        fh.write('<\html>\n') 

def _gen_tree_order_svgs_w_labs(svgs, labs, text_class = 'caption'):
    """
    Returns a single svg string from mupltipe svg strings and labels. 
    
    Parameters 
    ----------
    svgs : list of strs
        list of svgs strings
    name : str
        filename to output svg composite string
    dest : str
        destination folder
    text_class  : str
        CSS class name for labels useful for formating
    
    Returns 
    -------
    tree_order_svgs_w_labs : str

    """
    tree_order_svgs = list()
    for svg,lab in zip(svgs, labs):
        tree_order_svgs.append(f'<div><p class="{text_class}">{lab}</p></div>\n')
        tree_order_svgs.append(f"\t{svg}\n")
    tree_order_svgs_w_labs =  "\n".join(tree_order_svgs)
    return tree_order_svgs_w_labs

def _write_motif_tree_html( tree_order_svgs, 
                            name = "all_pal_svgs.html", 
                            template = 'cluster_tree_template.html',
                            dest = 'static',
                            write = True,
                            jinja_local = False):
    """
    Add tree_order_svgs string to a jinja2 html template 
    
    Parameters
    ----------

    tree_order_svgs : str
        string of multiple svgs separated by divs and label 
    name : str 
        name of the output file to be writen
    template : str
        template filename
    dest = 'static',
        output destination 
    write : bool
        If True write file.
    jinja_local : bool 
        If True, uses FileSystemLoader. If false, uses PackageLoader
    
    Returns
    ------
    tree_html : str

    """
    if jinja_local:
        # if you want to use a local directory called templates
        jinja_env = Environment(loader=FileSystemLoader('templates'))
    else:
        # if you want to use templates within the tcrdist2 package
        jinja_env = Environment(loader=PackageLoader('tcrdist', 'templates'))

    # the row template is for the motifs
    tree_tmp = jinja_env.get_template(template)
    # the page template is for the 
    tree_html = tree_tmp.render(tree_order_svgs  = tree_order_svgs)
    _write_html(html = tree_html, name =name  , dest = dest)

    return tree_html
    
def _modify_svg(svg, w ='500px'):
    """ Remove height and change width of an svg """
    svg = svg.replace('height="100%"', f'width="{w}"')
    svg = svg.replace('width="100%"', "")
    return svg  
