import numpy as np
import pandas as pd
import parasail
from tcrdist import pairwise
#from paths import path_to_matrices

class TCRrep:
    """
    Class for a TCRrep (T-Cell Receptor Repertoire) analysis.


    Attributes
    ----------
    cell_df : pandas.core.frame.DataFrame
        pandas.DataFrame containing data at the cell level
    clones_df: pandas.core.frame.DataFrame
        pandas.core.frame.DataFrame holding unique clones
    pwdist_df : pandas.core.frame.DataFrame
        pandas.DataFrame containing pairwise distances between unique unique_sequences
    index_cols : list
        list of strings, indicating columns in  for unique grouping
    meta_cols : list
        list of strings, indicating metadata columns in
    chains : list
        list of strings
    clones : None


    Methods
    -------
    deduplicate(self, index_cols)
        removes duplicate tcr-clones in
    compute_pairwise()
        computes pairwise distances on deduplicated data

    _validate_chains(self)
        raises ValueError is chains arg is mis-specified
    _validate_chain()
        raises ValueError if chain arg is mis-specified
    _validate_(self)
        raises TypeError if  is not pd.DataFrame
    _initialize_chain_specific_attributes(self)

    """

    def __init__(self, cell_df, chains=['alpha', 'beta']):
        self.cell_df = cell_df
        self.chains = chains
        self.pwdist_df = None
        self.clones_df = None
        self.index_cols = []
        self.meta_cols = None
        self.project_id = "<Your TCR Repertoire Project>"
        # VALIDATION OF INPUTS
        # check that chains are valid.
        self._validate_chains()
        # check that  is a pd.DataFrame
        self._validate_cell_df()
        # INITIALIZATION OF SPECIFIC ATTRIBUTES BASED ON SELECTED CHAINS
        self._initialize_chain_specific_attributes()


    def __repr__(self):
        return 'tcrdist.repertoire.TCRrep for {}\n with index_cols: {}\n'.format(self.project_id, self.index_cols)

    def __getitem__(self, position):
        # It should be decided whether get item should refer to the  or to the clone_df
        if self.clones_df is None:
            return self.cell_df.loc[position]
        if self.clones_df is not None:
            return self.clones_df.loc[position]

    def __len__(self):
        return self.cell_df.shape[0]

    def deduplicate(self):
        self.clones_df = _deduplicate(self.cell_df, self.index_cols)
        return self

    def compute_pairwise(self,
                         chain,
                         metric = "nw",
                         processes = 2,
                         user_function = None,
                         to_matrix = True,
                         **kwargs):

        # validate chain argument passed
        self._validate_chain(chain)
        # another option would be to loop through the a list of chains
        index_col_from_chain = {'alpha' : 'cdr3_a_aa',
                                'beta'  : 'cdr3_b_aa',
                                'gamma' : 'crd3_g_aa',
                                'delta' : 'cdr3_d_aa'}

        sequences = self.clones_df[index_col_from_chain[chain]]

        # Pull the default substitution matrix
        if chain == "alpha":
            smat = self.cdr3_a_aa_smat
        elif chain == "beta":
            smat = self.cdr3_b_aa_smat
        elif chain == 'gamma':
            smat = self.cdr3_g_aa_smat
        elif chain == "delta":
            smat = self.cdr3_d_aa_smat

        # If kwargs were passed use them, otherwise pass chain-sp. smat from above
        if ('matrix' in kwargs) or ("open" in kwargs):
            pw = _compute_pairwise(sequences = sequences,
                                   metric = metric,
                                   processes = processes,
                                   user_function = user_function,
                                   **kwargs)
        else:
            pw = _compute_pairwise(sequences = sequences,
                                   metric = metric,
                                   processes = processes,
                                   user_function = user_function,
                                   **{'matrix' : smat})


        if chain == "alpha":
            self.cdr3_a_aa_pw = pw
        elif chain == "beta":
            self.cdr3_b_aa_pw = pw
        elif chain == 'gamma':
            self.cdr3_g_aa_pw = pw
        elif chain == "delta":
            self.cdr3_d_aa_pw = pw


    def _validate_chains(self):
        """
        raise ValueError if invalid chains are passed to TCRrep __init__
        """
        check_chains_arg = ['alpha', 'beta', "gamma", "delta"]
        if len([c for c in self.chains if c not in check_chains_arg]) > 0:
            raise ValueError('TCRrep chains arg can be one or more of the'
                             'following {} case-sensitive'.format(check_chains_arg))

    def _validate_chain(self, chain):
        if chain not in ['alpha', 'beta', "gamma", "delta"]:
            raise ValueError('in compute_pairwise() chain must be one of the'
                             'following: "alpha", "beta", "gamma", "delta"' )

    def _validate_cell_df(self):
        """
        raise ValueError if  is not properly formatted.
        """
        if not isinstance(self.cell_df, pd.DataFrame):
            raise ValueError('TCRrep  argument must be pandas.DataFrame')
        # TODO: When know, validator should check column names and datatypes

    def _initialize_chain_specific_attributes(self):
        """
        Initialize pw object and default substitution matrix (smat) based on
        chains arguments.

        Naming of all objects have a standardized order
            region_chain_molecular_object
            (cdr3)_(a|b|d|g)_(aa|p)_(pw|smat|hmat)

        """
        if "alpha" in self.chains:
            self.cdr3_a_aa_pw = None
            self.cdr3_a_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_a_aa")
        if 'beta' in self.chains:
            self.cdr3_b_aa_pw = None
            self.cdr3_b_aa_smat= parasail.blosum62
            self.index_cols.append("cdr3_b_aa")

        if 'gamma' in self.chains:
            self.cdr3_g_aa_pw = None
            self.cdr3_g_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_g_aa")

        if 'delta' in self.chains:
            self.cdr3_d_aa_pw = None
            self.cdr3_d_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_d_aa")



    #@property
    #def clean_df(self):
    #        return self.[meta_cols + ['CDR3']]


def _deduplicate(cell_df, index_cols):
    clones = cell_df.groupby(index_cols)['count'].agg(np.sum).reset_index()
    return clones

def _compute_pairwise(sequences, metric = "nw", processes = 2, user_function = None, **kwargs):

    unique_seqs = pd.Series(sequences).unique()

    pw = pairwise.apply_pw_distance_metric_w_multiprocessing(
        sequences = sequences,
        metric = metric,
        user_function = user_function,
        processes= processes,
        **kwargs)

    pw = pairwise._pack_matrix(pw)
    pw_df = pd.DataFrame(pw, index = unique_seqs, columns = unique_seqs)
    pw_full = pw_df.loc[sequences, sequences]
    pw_full_np = pw_full.values
    return(pw_full_np)


# def _generate_parasail_hamming_aa_smat(filename = "inv_hamming.txt"):
#     """
#     Function that customizes a substitution matrix for use with parasail. It
#     customizes an existing blosum62 parasail substitution matrix to calculate
#     hamming distance.
#
#     Parameters
#     ----------
#     weight : int
#         integer to apply to all non diagonal entries
#
#     Returns
#     -------
#     x: parasail.bindings_v2.Matrix instance
#         returns parasail substitution matrix (approximating Hamming Distance)
#         with integer weight on all non-diagonal entries and zero on the diagonal.
#     """
#     return(parasail.Matrix( path_to_matrices + "/" + filename ))


# def _generate_parasail_hamming_aa_smat(weight = 1):
#     """
#     Parasail copy function as shown on github is not working so we have to import from
#     file.
#     Function that customizes a substitution matrix for use with parasail. It
#     customizes an existing blosum62 parasail substitution matrix to calculate
#     hamming distance.
#
#     Parameters
#     ----------
#     weight : int
#         integer to apply to all non diagonal entries
#
#     Returns
#     -------
#     x: parasail.bindings_v2.Matrix instance
#         returns parasail substitution matrix (approximating Hamming Distance)
#         with integer weight on all non-diagonal entries and zero on the diagonal.
#     """
#     x = parasail.blosum62.copy()
#     smat = np.ones(shape=(24,24), dtype= int)
#     np.fill_diagonal(smat, 0)
#     smat = np.multiply(smat, weight)
#     x.name = 'hamming_{}'.format(weight)
#     x.matrix = smat
#     return(x)
