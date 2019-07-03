import numpy as np
import pandas as pd
import parasail
from tcrdist import pairwise
from tcrdist.cdr3s_human import pb_cdrs
import warnings
#from paths import path_to_matrices

class TCRrep:
    """
    Class for a TCRrep (T-Cell Receptor Repertoire) analysis.


    Attributes
    ----------
    cell_df : pandas.core.frame.DataFrame
        pandas.DataFrame containing data at the cell level
    clone_df: pandas.core.frame.DataFrame
        pandas.core.frame.DataFrame holding unique clones
    pwdist_df : pandas.core.frame.DataFrame
        pandas.DataFrame containing pairwise distances between unique unique_sequences
    index_cols : list
        list of strings, indicating columns in  for unique grouping
    organism : string
        string either "human" or "mouse"
    meta_cols : list
        list of strings, indicating metadata columns in
    chains : list
        list of strings
    clones : None


    Methods
    -------
    infer_cdrs_from_v_gene()
        given v-gene name looks up amino acid sequences
    deduplicate(index_cols)
        removes duplicate tcr-clones in
    compute_pairwise_all()
        computes pairwise distances on deduplicated data for all

    Private Methods
    ---------------
    _validate_chains(self)
        raises ValueError is chains arg is mis-specified
    _validate_chain()
        raises ValueError if chain arg is mis-specified
    _validate_(self)
        raises TypeError if  is not pd.DataFrame
    _initialize_chain_specific_attributes(self, chain)
        creates chain specific attribute including setting default sub matrix
    _get_smat(self, chain, index_col)
        returns smat given chain (e.g. "alpha") and index_col (e.g. "cdr2_a_aa")
    _assign_pw_result(self, pw, chain, index_col)
        assigns pw distance given chain (e.g. "alpha") and index_col (e.g. "cdr2_a_aa")


    """

    def __init__(self, cell_df, chains=['alpha', 'beta'], organism = "human"):
        self.cell_df = cell_df
        self.chains = chains
        self.organism = organism
        self.pwdist_df = None
        self.clone_df = None
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
        return 'tcrdist.repertoire.TCRrep for {}\n with index_cols: {}\n with model organism: {}'.format(self.project_id, self.index_cols, self.organism)

    def __getitem__(self, position):
        # It should be decided whether get item should refer to the  or to the clone_df or it could be for iterating over pw dist matrices
        if self.clone_df is None:
            return self.cell_df.loc[position]
        if self.clone_df is not None:
            return self.clone_df.loc[position]

    def __len__(self):
        return self.cell_df.shape[0]

    def deduplicate(self):
        """
        With attribute self.index_col calls _deduplicate() and assigns
        result to attribute self.clone_df
        """
        self.clone_df = _deduplicate(self.cell_df, self.index_cols)
        return self

    def infer_cdrs_from_v_gene(self, chain):
        """
        Function which takes TCR V-gene names and infers the amino amino_acid
        sequence of key cdr1, cdr2, and pmhc regions.

        Parameters
    	----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'

        Assigns
    	-------
    	Assigns [cdr3|cdr2|cdr1|pmhc]_[a|b|d|g]_aa columns in self.cell_df

        Examples
    	-------
        testrep = TCRrep(cell_df = example_df, organism = "human", chains= ["alpha","beta"])
        testrep.infer_cdrs_from_v_gene(chain = "alpha")
        testrep.infer_cdrs_from_v_gene(chain = "beta")
        testrep.index_cols = testrep.index_cols + ['cdr1_a_aa','cdr2_a_aa', 'pmhc_a_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']

        Notes
    	-------
        Function which takes TCR V-gene names and infers the amino amino_acid
        sequence of key cdr1, cdr2, and pmhc region. (the pMHC-facing loop between
        CDR2 and CDR3 (IMGT alignment columns 81 - 86) based on lookup in
        dictionary: from tcrdist.cdr3s_human import pb_cdrs.

        IMGT Definitions of cdr1, cdr2, and pMHC-facing
        http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
        """
        f0 = lambda v : _map_gene_to_reference_seq(gene = v, cdr = 0)
        f1 = lambda v : _map_gene_to_reference_seq(gene = v, cdr = 1)
        f2 = lambda v : _map_gene_to_reference_seq(gene = v, cdr = 2)
        if chain is "alpha":
            self.cell_df['cdr1_a_aa'] = map(f0, self.cell_df.v_a_gene)
            self.cell_df['cdr2_a_aa'] = map(f1, self.cell_df.v_a_gene)
            self.cell_df['pmhc_a_aa'] = map(f2, self.cell_df.v_a_gene)
        if chain is "beta":
            self.cell_df['cdr1_b_aa'] = map(f0, self.cell_df.v_b_gene)
            self.cell_df['cdr2_b_aa'] = map(f1, self.cell_df.v_b_gene)
            self.cell_df['pmhc_b_aa'] = map(f2, self.cell_df.v_b_gene)
        if chain is "gamma":
            self.cell_df['cdr1_g_aa'] = map(f0, self.cell_df.v_g_gene)
            self.cell_df['cdr2_g_aa'] = map(f1, self.cell_df.v_g_gene)
            self.cell_df['pmhc_g_aa'] = map(f2, self.cell_df.v_g_gene)
        if chain is "delta":
            self.cell_df['cdr1_d_aa'] = map(f0, self.cell_df.v_d_gene)
            self.cell_df['cdr2_d_aa'] = map(f1, self.cell_df.v_d_gene)
            self.cell_df['pmhc_d_aa'] = map(f2, self.cell_df.v_d_gene)



    def compute_pairwise(self,
                         chain,
                         metric = "nw",
                         processes = 2,
                         user_function = None,
                         to_matrix = True,
                         **kwargs):
        """
        Early Function to be replaced with compute_pairwise_all
        """

        # validate chain argument passed
        self._validate_chain(chain)
        # another option would be to loop through the a list of chains
        index_col_from_chain = {'alpha' : 'cdr3_a_aa',
                                'beta'  : 'cdr3_b_aa',
                                'gamma' : 'crd3_g_aa',
                                'delta' : 'cdr3_d_aa'}

        sequences = self.clone_df[index_col_from_chain[chain]]

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

    def compute_pairwise_all(self,
                             chain,
                             compute_specific_region = None,
                             metric = "hamming",
                             processes = 2,
                             user_function = None,
                             to_matrix = True,
                             **kwargs):
        """
        Function that computes pairwise distances for all regions on a given
        chain or for a specific region.

        Parameters
    	----------
        chain: string
            'alpha', 'beta', 'gamma', or 'delta'
        compute_specific_region : string
            string (e.g. "cdr2_a_aa") to over-ride function behavior and compute
            only a single region
        metric : string
            'nw', 'hamming', or 'custom'
        processes: int
            int for number of available cpu for multiprocessing (to see available
            try multiprocessing.cpu_count())
        user_function: function
            function for a custom distance metric on two strings (This is
            an advanced option, so don't use this unless you are absolutely
            sure what you are doing; metric arg must be set to 'custom').
        to_matrix: boolean
            True will return pairwise distance as result as a 2D ndarray

    	Assigns
    	-------
    	self.[cdr3|cdr2|cdr1|pmhc]_[a|b|d|g]_aa_pw objects

        Example
    	-------
        testrep = TCRrep(cell_df = example_df, organism = "human", chains= ["alpha","beta"])
        testrep.infer_cdrs_from_v_gene(chain = "alpha")
        testrep.infer_cdrs_from_v_gene(chain = "beta")
        testrep.index_cols = testrep.index_cols + ['cdr1_a_aa','cdr2_a_aa','pmhc_a_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']
        testrep.deduplicate()
        testrep.compute_pairwise_all(chain = "alpha", metric= "hamming")
        testrep.compute_pairwise_all(chain = "beta", metric= "hamming")

        # alternatively, compute each region one by one
        testrep.compute_pairwise_all(chain = "beta", compute_specific_region="cdr1_b_aa")
        testrep.compute_pairwise_all(chain = "alpha", compute_specific_region="cdr2_a_aa")

        """

        # validate chain argument passed
        self._validate_chain(chain)
        # If compute_specific_region is None, then the behavior is to loop through the a list regions.
        if compute_specific_region is None:
            index_col_from_chain = {'alpha' : ['cdr3_a_aa','cdr2_a_aa',
                                               'cdr1_a_aa','pmhc_a_aa'],
                                    'beta'  : ['cdr3_b_aa','cdr2_b_aa',
                                               'cdr1_b_aa','pmhc_b_aa'],
                                    'gamma' : ['cdr3_g_aa','cdr2_g_aa',
                                               'cdr1_g_aa','pmhc_g_aa'],
                                    'delta' : ['cdr3_d_aa','cdr2_d_aa',
                                               'cdr1_d_aa','pmhc_d_aa']}
        # Alternative behavior: is to loop over a single chain and region.
        else:
            index_col_from_chain = {}
            index_col_from_chain[chain] = [compute_specific_region]


        for index_col in index_col_from_chain[chain]:
            try:
                sequences = self.clone_df[index_col]
            except KeyError:
                warnings.warn("{} not found, no distances computed for {}".format(index_col, index_col))
                continue



            # COMPUTE PAIRWISE
            # If kwargs were passed use them, otherwise pass chain-sp. smat from above
            if ('matrix' in kwargs) or ("open" in kwargs):
                pw = _compute_pairwise(sequences = sequences,
                                       metric = metric,
                                       processes = processes,
                                       user_function = user_function,
                                       **kwargs)
            else:
                # Pull the default substitution matrix from object attributes
                smat = self._get_smat(chain = chain, index_col = index_col)
                pw = _compute_pairwise(sequences = sequences,
                                       metric = metric,
                                       processes = processes,
                                       user_function = user_function,
                                       **{'matrix' : smat})

            # ASSIGN RESULT
            self._assign_pw_result(pw = pw, chain=chain, index_col=index_col)

    def compute_paired_tcrdist(self,
                               chains = ['alpha', 'beta'],
                               replacement_weights = {}):
        """
        Computes tcrdistance metric combining distances metrics across multiple
        TCR regions.

        Parameters
        ----------
        chains: list
            list of strings containing some combination of 'alpha', 'beta',
            'gamma', and 'delta'
        replacement_weights: dictionary
            optional dictionary of the form {'cdr1_a_aa_pw':1, 'cdr2_a_aa_pw':1}
            used to place greater weight on certain TCR regions. The default
            is a weight of 1.

        Returns
        -------
        r: dictionary
            dictionary with a 2D tcrdist np.array and dictionary of regions and
            relative weights
        """
        [self._validate_chain(c) for c in chains]
        weights = {'cdr1_a_aa_pw':1,
                   'cdr2_a_aa_pw':1,
                   'cdr3_a_aa_pw':1,
                   'pmhc_a_aa_pw':1,
                   'cdr1_b_aa_pw':1,
                   'cdr2_b_aa_pw':1,
                   'cdr3_b_aa_pw':1,
                   'pmhc_b_aa_pw':1,
                   'cdr1_g_aa_pw':1,
                   'cdr2_g_aa_pw':1,
                   'cdr3_g_aa_pw':1,
                   'pmhc_g_aa_pw':1,
                   'cdr1_d_aa_pw':1,
                   'cdr2_d_aa_pw':1,
                   'cdr3_d_aa_pw':1,
                   'pmhc_d_aa_pw':1}

        for k in replacement_weights:
            weights[k] = replacement_weights[k]

        alpha_keys = [k for k in weights.keys() if k.endswith("a_aa_pw")]
        beta_keys  = [k for k in weights.keys() if k.endswith("b_aa_pw")]
        gamma_keys = [k for k in weights.keys() if k.endswith("g_aa_pw")]
        delta_keys = [k for k in weights.keys() if k.endswith("d_aa_pw")]

        full_keys = []
        if 'alpha' in chains:
            full_keys = full_keys + alpha_keys
        if 'beta' in chains:
            full_keys = full_keys + beta_keys
        if 'gamma' in chains:
            full_keys = full_keys + gamma_keys
        if 'delta' in chains:
            full_keys = full_keys + delta_keys

        tcrdist = np.zeros(self.__dict__[full_keys[0]].shape)
        for k in full_keys:
            tcrdist = self.__dict__[k]*weights[k] + tcrdist
        self.paired_tcrdist = tcrdist
        self.paired_tcrdist_weights = {k:weights[k] for k in full_keys}
        r = {'paired_tcrdist' : tcrdist,
                'paired_tcrdist_weights' : {k:weights[k] for k in full_keys}}
        return(r)

    def _validate_chains(self):
        """
        raise ValueError if invalid chains are passed to TCRrep __init__
        """
        check_chains_arg = ['alpha', 'beta', "gamma", "delta"]
        if len([c for c in self.chains if c not in check_chains_arg]) > 0:
            raise ValueError('TCRrep chains arg can be one or more of the '
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
            self.cdr3_a_aa_smat = parasail.blosum62
            self.cdr2_a_aa_smat = parasail.blosum62
            self.cdr1_a_aa_smat = parasail.blosum62
            self.pmhc_a_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_a_aa")

        if 'beta' in self.chains:
            self.cdr3_b_aa_smat = parasail.blosum62
            self.cdr2_b_aa_smat = parasail.blosum62
            self.cdr1_b_aa_smat = parasail.blosum62
            self.pmhc_b_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_b_aa")

        if 'gamma' in self.chains:
            self.cdr3_g_aa_smat = parasail.blosum62
            self.cdr2_g_aa_smat = parasail.blosum62
            self.cdr1_g_aa_smat = parasail.blosum62
            self.pmhc_g_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_g_aa")

        if 'delta' in self.chains:
            self.cdr3_d_aa_smat = parasail.blosum62
            self.cdr2_d_aa_smat = parasail.blosum62
            self.cdr1_d_aa_smat = parasail.blosum62
            self.pmhc_d_aa_smat = parasail.blosum62
            self.index_cols.append("cdr3_d_aa")



    def _get_smat(self, chain, index_col):
        """
        Gets the correct substitution matrix (smat) based on chain and column

        Parameters
        ----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        index_col: string
            [cdr3|cdr2|cdr1|pmhc]_[a|b|g|d]_aa_pw
        """
        self._validate_chain(chain = chain)

        if chain == "alpha":
            if index_col.startswith("cdr3_a"):
                smat = self.cdr3_a_aa_smat
            elif index_col.startswith("cdr2_a"):
                smat = self.cdr2_a_aa_smat
            elif index_col.startswith("cdr1_a"):
                smat = self.cdr1_a_aa_smat
            elif index_col.startswith("pmhc_a"):
                smat = self.pmhc_a_aa_smat
            else:
                smat = parasail.blosum62
                warnings.warn("Using default parasail.blosum62 because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))
        if chain == "beta":
            if index_col.startswith("cdr3_b"):
                smat = self.cdr3_b_aa_smat
            elif index_col.startswith("cdr2_b"):
                smat = self.cdr2_b_aa_smat
            elif index_col.startswith("cdr1_b"):
                smat = self.cdr1_b_aa_smat
            elif index_col.startswith("pmhc_b"):
                smat = self.pmhc_b_aa_smat
            else:
                smat = parasail.blosum62
                warnings.warn("Using default parasail.blosum62 because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))
        if chain == "gamma":
            if index_col.startswith("cdr3_g"):
                smat = self.cdr3_g_aa_smat
            elif index_col.startswith("cdr2_g"):
                smat = self.cdr2_g_aa_smat
            elif index_col.startswith("cdr1_g"):
                smat = self.cdr1_g_aa_smat
            elif index_col.startswith("pmhc_g"):
                smat = self.pmhc_g_aa_smat
            else:
                smat = parasail.blosum62
                warnings.warn("Using default parasail.blosum62 because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))
        if chain == "delta":
            if index_col.startswith("cdr3_d"):
                smat = self.cdr3_d_aa_smat
            elif index_col.startswith("cdr2_d"):
                smat = self.cdr2_d_aa_smat
            elif index_col.startswith("cdr1_d"):
                smat = self.cdr1_d_aa_smat
            elif index_col.startswith("pmhc_d"):
                smat = self.pmhc_d_aa_smat
            else:
                smat = parasail.blosum62
                warnings.warn("Using default parasail.blosum62 because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))

        return(smat)


    def _assign_pw_result(self, pw, chain, index_col):
        """
        Assigns pairwise result to TCRrep attribute based on chain and index_col

        Parameters
        ----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        index_col: string
            [cdr3|cdr2|cdr1|pmhc]_[a|b|g|d]_aa_pw

        """
        self._validate_chain(chain = chain)

        if chain == "alpha":
            if index_col.startswith("cdr3_a"):
                self.cdr3_a_aa_pw = pw
            elif index_col.startswith("cdr2_a"):
                self.cdr2_a_aa_pw = pw
            elif index_col.startswith("cdr1_a"):
                self.cdr1_a_aa_pw = pw
            elif index_col.startswith("pmhc_a"):
                self.pmhc_a_aa_pw = pw
            else:
                warnings.warn("No assignment for {} because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))

        elif chain == "beta":
            if index_col.startswith("cdr3_b"):
                self.cdr3_b_aa_pw = pw
            elif index_col.startswith("cdr2_b"):
                self.cdr2_b_aa_pw = pw
            elif index_col.startswith("cdr1_b"):
                self.cdr1_b_aa_pw = pw
            elif index_col.startswith("pmhc_b"):
                self.pmhc_b_aa_pw = pw
            else:
                warnings.warn("No assignment for {} because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))

        elif chain == 'gamma':
            if index_col.startswith("cdr3_g"):
                self.cdr3_g_aa_pw = pw
            elif index_col.startswith("cdr2_g"):
                self.cdr2_g_aa_pw = pw
            elif index_col.startswith("cdr1_g"):
                self.cdr1_g_aa_pw = pw
            elif index_col.startswith("pmhc_g"):
                self.pmhc_g_aa_pw = pw
            else:
                warnings.warn("No assignment for {} because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))

        elif chain == "delta":
            if index_col.startswith("cdr3_d"):
                self.cdr3_d_aa_pw = pw
            elif index_col.startswith("cdr2_d"):
                self.cdr2_d_aa_pw = pw
            elif index_col.startswith("cdr1_d"):
                self.cdr1_d_aa_pw = pw
            elif index_col.startswith("pmhc_d"):
                self.pmhc_d_aa_pw = pw
            else:
                warnings.warn("No assignment for {} because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))




def _map_gene_to_reference_seq(organism = "human",
                               gene= 'TRAV1-1*02',
                               cdr = 1,
                               ref = pb_cdrs):
    """
    Function that takes TCR V-gene


    Parameters
    ---------
    organism : string
        string must be "human" or "mouse"

    gene : string
        string specifying gene (e.g 'TRAV1-1*02)

    Returns
    -------
    aa_string : string
        amino acid string or None if gene not in ref

    """
    try:
        aa_string = ref[organism][gene][cdr][0]
    except KeyError:
        aa_string = None
    return aa_string

def _deduplicate(cell_df, index_cols):
    clones = cell_df.groupby(index_cols)['count'].agg(np.sum).reset_index()
    return clones

def _compute_pairwise(sequences, metric = "nw", processes = 2, user_function = None, **kwargs):

    unique_seqs = pd.Series(sequences).unique()

    pw = pairwise.apply_pw_distance_metric_w_multiprocessing(
        sequences = unique_seqs, #! BUG FIX (sequences changed to unique_seqs)
        metric = metric,
        user_function = user_function,
        processes= processes,
        **kwargs)

    pw = pairwise._pack_matrix(pw)
    pw_df = pd.DataFrame(pw, index = unique_seqs, columns = unique_seqs)
    pw_full = pw_df.loc[sequences, sequences]
    pw_full_np = pw_full.values
    return(pw_full_np)
