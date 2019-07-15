"""
mappers
=======

The mappers module contain dictionaries and methods for taking common output and
reformatting it to match tcrdist2 input requirements.

The top of the file contains ordered dictionaries mapping keys

More complex mappings are encapsulated in methods to be called on pd.DataFrames.

"""
import ast
import pandas as pd
from collections import OrderedDict
from functools import reduce

tcrdist_to_tcrdist2_mapping = OrderedDict([('id', 'id'),
                                           ('epitope', 'epitope'),
                                           ('subject', 'subject'),
                                           ('cdr3a', 'cdr3_a_aa'),
                                           ('cdr3b', 'cdr3_b_aa'),
                                           ('ja_gene', 'j_a_gene'),
                                           ('va_gene', 'v_a_gene'),
                                           ('jb_gene', 'j_b_gene'),
                                           ('vb_gene', 'v_b_gene')])


vdjdb_to_tcrdist2_mapping_TRA = OrderedDict([('complex.id', 'complex_id'),
                                        ('Gene', 'gene'),
                                        ('CDR3', 'cdr3_a_aa'),
                                        ('V', 'v_a_gene'),
                                        ('J', 'j_a_gene'),
                                        ('Species', 'organism'),
                                        ('MHC A', 'mhc_a_a'),
                                        ('MHC B', 'mhc_a_b'),
                                        ('MHC class', 'mhc_a_class'),
                                        ('Epitope', 'epitope'),
                                        ('Epitope gene', 'epitope_gene'),
                                        ('Epitope species', 'epitope_species'),
                                        ('Reference', 'reference'),
                                        ('Method', 'method'),
                                        ('Meta', 'meta'),
                                        ('CDR3fix', 'cdr3fix'),
                                        ('Score', 'score')])

vdjdb_to_tcrdist2_mapping_TRB = OrderedDict([('complex.id', 'complex_id'),
                                        ('Gene', 'gene'),
                                        ('CDR3', 'cdr3_b_aa'),
                                        ('V', 'v_b_gene'),
                                        ('J', 'j_b_gene'),
                                        ('Species', 'organism'),
                                        ('MHC A', 'mhc_b_a'),
                                        ('MHC B', 'mhc_b_b'),
                                        ('MHC class', 'mhc_b_class'),
                                        ('Epitope', 'epitope'),
                                        ('Epitope gene', 'epitope_gene'),
                                        ('Epitope species', 'epitope_species'),
                                        ('Reference', 'reference'),
                                        ('Method', 'method'),
                                        ('Meta', 'meta'),
                                        ('CDR3fix', 'cdr3fix'),
                                        ('Score', 'score')])




def vdjdb_to_tcrdist2(pd_df):
    """
    reformat .tsv downloaded from vdjdb (July 2019) to input format of tcrdist2 input format

    Parameters
    ----------
    pd_df : pandas.core.DataFrame likley derived from pd.read_csv() on .tsv from
    vdjdb.org

    Returns
    -------
    tcrdist2_formatted_pd_df

    Notes
    -----
    Code is written for clarity rather than length or computatinal efficiency
    since parsing demands may change.

    Note that in the unpack phase in the fist block of code a dictionary is made
    to hold all features in the vdjdb output. Much of the information from
    paired sequences is redundant (i.e. epitope). If memory demands increase
    for larger number of rows, consider modifying to unpack less of the output.

    Requires: pandas and ast (abstract syntax trees)

    Raises
    ------
    AssertionError if 'Gene' field is not TRA and TRB

    TODO
    ----
    Handle case of missing alpha or beta seq


    """
    d = {}
    # unpack vdjdb panda.DataFrame to dictionary
    for index, row in pd_df.iterrows():

        # gene will be either TRA or
        gene = row['Gene']
        assert(gene in ["TRA", "TRB"]), "Unexpected Input in vdjdb 'Gene' field not TRA or TRB. Check input."

        # complex_id is shared for paired chains
        complex_id = row['complex.id']

        # note: ast (abstract syntax trees) allows converstion of string to dictionary)
        # packed within a pd.DF cell
        meta_dict =  ast.literal_eval(row['Meta'])
        method_dict = ast.literal_eval(row['Method'])

        if gene == "TRA":
            row = row.rename(vdjdb_to_tcrdist2_mapping_TRA)
        elif gene == "TRB":
            row = row.rename(vdjdb_to_tcrdist2_mapping_TRB)

        d.setdefault(complex_id, {})[gene] = row.to_dict()
        d[complex_id][gene].update(meta_dict)
        d[complex_id][gene].update(method_dict)


    # output select fields to a list of dictionaries (l_out)
    complex_ids = sorted(d.keys())
    l_out = []
    for complex_id in complex_ids:
        try:
            id         = d[complex_id]["TRA"]['complex_id']
            cell_type  = d[complex_id]["TRA"]['cell.subset']
            organism   = d[complex_id]["TRA"]['organism']
            epitope_aa = d[complex_id]["TRA"]['epitope']
            epitope    = d[complex_id]["TRA"]['epitope.id']
            subject    = d[complex_id]["TRA"]['subject.id']

            mhc_a_a    = d[complex_id]["TRA"]['mhc_a_a']
            mhc_a_b    = d[complex_id]["TRA"]['mhc_a_b']

            mhc_b_a    = d[complex_id]["TRB"]['mhc_b_a']
            mhc_b_b    = d[complex_id]["TRB"]['mhc_b_b']

            cdr3_a_aa  = d[complex_id]["TRA"]['cdr3_a_aa']
            v_a_gene   = d[complex_id]["TRA"]['v_a_gene']
            j_a_gene   = d[complex_id]["TRA"]['j_a_gene']

            cdr3_b_aa  = d[complex_id]["TRB"]['cdr3_b_aa']
            v_b_gene   = d[complex_id]["TRB"]['v_b_gene']
            j_b_gene   = d[complex_id]["TRB"]['j_b_gene']

            frequency  = d[complex_id]["TRA"]["frequency"]
            try:
                count      = int(frequency.split("/")[0])
            except ValueError:
                count      = 1
            l_out.append(
            {'id'        : id,
            'cell_type'  : cell_type,
            'organism'   : organism,
            'epitope_aa' : epitope_aa,
            'epitope'    : epitope,
            'subject'    : subject,
            'mhc_a_a'    : mhc_a_a,
            'mhc_a_b'    : mhc_a_b,
            'mhc_b_a'    : mhc_b_a,
            'mhc_b_b'    : mhc_b_b,
            'cdr3_a_aa'  : cdr3_a_aa,
            'v_a_gene'   : v_a_gene,
            'j_a_gene'   : j_a_gene,
            'cdr3_b_aa'  : cdr3_b_aa,
            'v_b_gene'   : v_b_gene,
            'j_b_gene'   : j_b_gene,
            'frequency'  : frequency,
            'count'      : count})
        except KeyError:
            pass


    # convert list of dictionaries to pandas DataFrame
    tcrdist2_formatted_pd_df = pd.DataFrame.from_dict(l_out)
    return(tcrdist2_formatted_pd_df)



def multistudy_vdjdb_to_tcrdist2(pd_df):
    """

    map vdjdb data from multiple studies to one tcrdist2 formatted pd.dataframe

    Parameters
    ----------
    pd_df : Pandas DataFrame
        from vdjdb.org with headers as of July 10, 2019

    Returns
    -------

    tcr_formatted_dfs_merged : Pandas DataFrame
        with all paired data from input, single combined data frame tcrdist2 mapped

    """
    # 2 complex.id which links paired sequences is Reference (study specific)
    studies = pd_df['Reference'].unique()

    # break full df into sub dfs: one per study
    dfs_split_by_study = {study: pd.DataFrame for study in studies}
    for study in list(dfs_split_by_study.keys()):
        dfs_split_by_study[study] = pd_df[:][pd_df.Reference == study].reset_index()

    # tcrdist format each study specific dataframe
    tcr_formatted_dfs_split_by_study = {study: None for study in studies}
    for study in list(dfs_split_by_study.keys()):
        tcr_formatted_dfs_split_by_study[study] = vdjdb_to_tcrdist2(pd_df = dfs_split_by_study[study])

    # create list of studies with paired data
    studies_with_AB_data = [k for k in list(tcr_formatted_dfs_split_by_study.keys()) if \
                            tcr_formatted_dfs_split_by_study[k].shape[0] is not 0]
    studies_without_AB_data = [k for k in list(tcr_formatted_dfs_split_by_study.keys()) if \
                               tcr_formatted_dfs_split_by_study[k].shape[0] is not 0]

    # subset dictionary of dataframes to those with paired data
    ab = { k: tcr_formatted_dfs_split_by_study[k] for k in studies_with_AB_data }

    # count rows in each dataset for human interest
    ab_rows = [x[1].shape[0] for x in list(ab.items())]

    # convert dictionary of dataframes into a list and only keep the dataframe ignoring the study name.
    def append_study_name(name, df):
        df["reference"] = name
        return(df)

    ab_dfs = [append_study_name(x[0], x[1]) for x in list(ab.items())]

    # combine the list of dataframes into a single dataframe
    tcr_formatted_dfs_merged = reduce(lambda  x, y: pd.concat([x, y]), ab_dfs)

    return tcr_formatted_dfs_merged
