import numpy as np
import pandas as pd
import tarfile
import sys
import os
import scipy.spatial
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import collections
import json
import warnings
import pickle
import multiprocessing
import parasail
import pwseqdist
from zipdist.zip2 import Zipdist2 

from . import repertoire_db
from . import pgen
from . import mappers
from . import pairwise

# includes tools for use with explore.py
#from paths import path_to_matrices
#This replaces: from tcrdist.cdr3s_human import pb_cdrs
pb_cdrs = repertoire_db.generate_pbr_cdr()

class TCRrep:
    """
    Class for managing a T-Cell Receptor Repertoire (TCRrep) analysis. Produce
    a distance measure based on comparisons from multiple T-Cell receptor
    complementarity-determining regions (CDRs)


    Attributes
    ----------
    cell_df : pandas.core.frame.DataFrame
        input data at the level of individual cell level
    clone_df : pandas.core.frame.DataFrame
        deduplicated data frame at the level of unique clones
    index_cols : list
        list of strings, indicating columns to group cells to clones
    organism : string
        either "human" or "mouse"
    meta_cols : list
        list of strings, indicating metadata columns (e.g. hla_type)
    chains : list
        list of strings containing one or more of 'alpha', 'beta', 'gamma' or 'delta'
    stored_tcrdist : list
        list containing all previously generated outputs of
        `TCRrep.compute_paired_tcrdist`
    paired_tcrdist : ndarray
        most recent output of :py:meth:`tcrdist.repertoire.TCRrep.compute_paired_tcrdist`
    paired_tcrdist_weights : dictionary
        CDR weights used to generate the most recent output of
        TCRrep.compute_paired_tcrdist`
    all_genes : dictionary
        dictionary of reference TCRs

    Methods
    -------
    TCRrep.infer_cdrs_from_v_gene()
        infer CDR amino acid sequences from v-gene specified
    deduplicate()
        remove duplicate clones by grouping
    compute_pairwise_all()
        compute pairwise distances on deduplicated data for all regions in
        a chain. Alternatively can compute distance between a
    compute_paired_tcrdist()
        calculate weighted pairwise distance across all CDRs
    generate_ref_genes_from_db()
        generates all_genes attribute a dictionary of reference TCRs


    """
    def __init__(self,
                 cell_df,
                 chains=['alpha', 'beta'],
                 organism = "human",
                 db_file = "alphabeta_db.tsv"):
        self.db_file = db_file
        self.cell_df = cell_df
        self.chains = chains
        self.organism = organism
        self.pwdist_df = None
        self.clone_df = None
        self.index_cols = []
        self.stored_tcrdist = []
        self.paired_tcrdist = None
        self.paired_tcrdist_weights = None
        self.meta_cols = None
        self.project_id = "<Your TCR Repertoire Project>"
        self.all_genes = None
        self.imgt_aligned_status = None
        
        # VALIDATION OF INPUTS
        # check that chains are valid.
        self._validate_organism()
        self._validate_chains()
        # check that  is a pd.DataFrame
        self._validate_cell_df()

        # INIT OF SPECIFIC ATTRIBUTES BASED ON SELECTED CHAINS
        self._initialize_chain_specific_attributes()
        # INIT the REFERENCE DB see repertoire_db.py
        self.generate_ref_genes_from_db(db_file)



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

    def generate_ref_genes_from_db(self, db_file = "alphabeta_db.tsv"):
        """

        Responsible for generating the all_genes attribute containing all
        the reference TCR data.

        Parameters
        ----------
        db_file : string

        Returns an ordered dictionary of reference sequences

        """
        self.all_genes = repertoire_db.RefGeneSet(db_file).all_genes

    def _map_gene_to_reference_seq2(self,
                                    organism,
                                    gene,
                                    cdr,
                                    attr ='cdrs_no_gaps'):
        """
        internal function that looks up the cdr sequence (gapped or ungapped)
        from the self.all_genes library

        Parameter
        ---------

        organism : string
            mouse or human
        gene : string
            specifies the TCR gene such as 'TRAV1*01'
        cdr : int
            0 - CDR1, 1-CDR2 and 2 - CDR2.5
        attr : string
            'cdrs_no_gaps' or 'cdrs_aligned' with gaps from IMGT
        """
        try:
            aa_string = self.all_genes[organism][gene].__dict__[attr][cdr]
        except KeyError:
            aa_string = None
            warnings.warn("{} gene was not recognized in reference db no cdr seq could be inferred".format(gene))
        return(aa_string)


    def deduplicate(self):
        """
        With attribute self.index_col calls _deduplicate() and assigns
        result to attribute self.clone_df
        """
        self.clone_df = _deduplicate(self.cell_df, self.index_cols)
        
        # check if any clones were lost due to missing information
        if np.sum(self.cell_df['count']) != np.sum(self.clone_df['count']):
            n_cells_lost = np.sum(self.cell_df['count']) - np.sum(self.clone_df['count'])
            n_cell = np.sum(self.cell_df['count'])
            warnings.warn(f"Not all cells/sequences could be grouped into clones. {n_cells_lost} of {n_cell} were not captured. This occurs when any of the values in the index columns are null or missing for a given sequence. To see entries with missing values use: tcrdist.repertoire.TCRrep.show_incomplete()\n")
        
        # if no clone id column provided thetrn create one as a sequence of numbers
        if "clone_id" not in self.clone_df:
            N = self.clone_df.shape[0]
            self.clone_df['clone_id'] = range(1, N + 1 ,1)
        
        return self

    def show_incomplete(self):      
        ind = self.cell_df[self.index_cols].isnull().any(axis = 1)   
        incomplete_clones = self.cell_df.loc[ind,self.index_cols].copy()
        return incomplete_clones  

    # def tcr_motif_clones_df(self):
    #     """
    #     Use this function to create a clones_df input appropriate to TCRMotif.
    #
    #     It make use of a mapper to ensure proper columns and column names
    #
    #     Example
    #     -------
    #     TCRMotif(clones_df = TCRRep.tcr_motif_clones_df())
    #     """
    #     return _map_clone_df_to_TCRMotif_clone_df(self.clone_df)

    def tcr_motif_clones_df(self):
        """
        Use this function to create a clones_df input appropriate to TCRMotif.

        It make use of a mapper to ensure proper columns and column names

        Example
        -------
        TCRMotif(clones_df = TCRrep.tcr_motif_clones_df())
        """
        return mappers.generic_pandas_mapper(self.clone_df,
                                             mappers.TCRrep_clone_df_to_TCRMotif_clone_df)


    def infer_cdrs_from_v_gene(self, chain, imgt_aligned = False):
        """
        Function taking TCR v-gene name to infer the amino amino_acid
        sequence of cdr1, cdr2, and pmhc loop regions.

        Parameters
    	----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        imgt_aligned : boolean
            if True cdr1, cdr2, cdr2.5 will be returned with gaps
            and by definition will be the same length. MSH.......ET


        Returns
    	-------
        self.cell_df : pandas.core.frame.DataFrame
    	   Assigns [cdr3|cdr2|cdr1|pmhc]_[a|b|d|g]_aa columns in self.cell_df

        Examples
    	--------
        >>> testrep = TCRrep(cell_df = example_df, organism = "human", chains= ["alpha","beta"])
        >>> testrep.infer_cdrs_from_v_gene(chain = "alpha")
        >>> testrep.infer_cdrs_from_v_gene(chain = "beta")
        >>> testrep.index_cols = testrep.index_cols + ['cdr1_a_aa','cdr2_a_aa', 'pmhc_a_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']

        Notes
    	-----
        This function takes the V-gene names and infers the amino amino_acid
        sequence of the cdr1, cdr2, and pmhc region (pmhc refers to the
        pMHC-facing loop between CDR2 and CDR3 (IMGT alignment columns 81 - 86.
        These sequences are based up on lookup from the dictionary here:

        originally: from tcrdist.cdr3s_human import pb_cdrs

        now:

        self.generate_ref_genes_from_db(db_file)

        imgt_aligned : boolean
            if True cdr1, cdr2, cdr2.5 will be returned with gaps
            and by definition will be the same length.
            MSH.......ET
            FNH.......DT
            LGH.......NA

        References
        ----------

        IMGT definitions of cdr1, cdr2, and pMHC-facing can be found here
        http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
        """

        if not imgt_aligned:
            self.imgt_aligned_status = False
            f0 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 0,
                                                             organism = self.organism,
                                                             attr ='cdrs_no_gaps')
            f1 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 1,
                                                             organism = self.organism,
                                                             attr ='cdrs_no_gaps')
            f2 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 2,
                                                             organism = self.organism,
                                                             attr ='cdrs_no_gaps')
        else:
            self.imgt_aligned_status = True
            f0 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 0,
                                                             organism = self.organism,
                                                             attr ='cdrs')
            f1 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 1,
                                                             organism = self.organism,
                                                             attr ='cdrs')
            f2 = lambda v : self._map_gene_to_reference_seq2(gene = v,
                                                             cdr = 2,
                                                             organism = self.organism,
                                                             attr ='cdrs')
        if chain is "alpha":
            self.cell_df['cdr1_a_aa'] = list(map(f0, self.cell_df.v_a_gene))
            self.cell_df['cdr2_a_aa'] = list(map(f1, self.cell_df.v_a_gene))
            self.cell_df['pmhc_a_aa'] = list(map(f2, self.cell_df.v_a_gene))
        if chain is "beta":
            self.cell_df['cdr1_b_aa'] = list(map(f0, self.cell_df.v_b_gene))
            self.cell_df['cdr2_b_aa'] = list(map(f1, self.cell_df.v_b_gene))
            self.cell_df['pmhc_b_aa'] = list(map(f2, self.cell_df.v_b_gene))
        if chain is "gamma":
            self.cell_df['cdr1_g_aa'] = list(map(f0, self.cell_df.v_g_gene))
            self.cell_df['cdr2_g_aa'] = list(map(f1, self.cell_df.v_g_gene))
            self.cell_df['pmhc_g_aa'] = list(map(f2, self.cell_df.v_g_gene))
        if chain is "delta":
            self.cell_df['cdr1_d_aa'] = list(map(f0, self.cell_df.v_d_gene))
            self.cell_df['cdr2_d_aa'] = list(map(f1, self.cell_df.v_d_gene))
            self.cell_df['pmhc_d_aa'] = list(map(f2, self.cell_df.v_d_gene))


    def infer_olga_aa_cdr3_pgens(self,
                                 chain,
                                 cdr3_only = False,
                                 chain_folder = None,
                                 recomb_type = None):
        """
        Infer the probability of generation using the Olga Code base
        (Sethna et al. 2018) updated to python 3 for use with tcrdist.

        Parameters
        ----------
        chain : string
            'alpha', 'beta' (TODO: create default models for 'gamma' and 'delta')
        cdr3_only : boolean
            (optional) if True, the amino acid cdr3 probability of generation statistic
            will be calculated without using the V or J gene usage statistics
        chain_folder : string
            (optional) specifies the OLGA default model folder containing a
            generative model. When None (which is recommended), the default
            folder is chosen based on the chain argument.
        recomb_type : string
            (optional) 'VDJ' or 'VJ' specifying the OLGA recombination model.
            When None (which is recommended), the default folder is chosen based
            on the chain argument.

        Returns
        -------
        olga_pgens : pd.Series
            containing the probability of generation, this output is also assigned
            to clone_df.cdr3_[a|b|g|d]_aa_pgen

        Notes
        -----
        tcrdist2 authors UPDATED THE FOLLOWING CODE TO PYTHON 3
        USING COMMIT e825c333f0f9a4eb02132e0bcf86f0dca9123114 (Jan 18, 2019)

        ORIGINAL OLGA CODE CAN BE FOUND AT:
        https://github.com/zsethna/OLGA

        """

        assert(isinstance(self.clone_df, pd.DataFrame)), "this function requires a valid TCRrep.clone_df has been instantiated"

        # The Nested If Statements assigns cdr3s, v_genes, j_genes based on chain, organism and other optional args
        if chain == "alpha":
            if (chain_folder is None):
                if self.organism is 'human':
                    chain_folder = "human_T_alpha"
                elif self.organism is 'mouse':
                    raise ValueError("SORRY: OLGA default files do not yet support mouse alpha TCRs")
                    chain_folder = "mouse_T_alpha"
            if (recomb_type is None):
                recomb_type = "VJ"
            cdr3s = self.clone_df.cdr3_a_aa

            if not cdr3_only:
                v_genes = self.clone_df.v_a_gene
                j_genes = self.clone_df.j_a_gene
            else:
                v_genes = None
                j_genes = None

        if chain == "beta":
            if (chain_folder is None):
                if self.organism is 'human':
                    chain_folder = "human_T_beta"
                elif self.organism is 'mouse':
                    chain_folder = "mouse_T_beta"
            if (recomb_type is None):
                recomb_type = "VDJ"
            cdr3s = self.clone_df.cdr3_b_aa

            if not cdr3_only:
                v_genes = self.clone_df.v_b_gene
                j_genes = self.clone_df.j_b_gene
            else:
                v_genes = None
                j_genes = None

        if chain  ==  "gamma":
            raise ValueError("SORRY: OLGA default files do not yet support gamma TCRs")
            if (chain_folder is None):
                if self.organism is 'human':
                    chain_folder = "human_T_gamma"
                elif self.organism is 'mouse':
                    chain_folder = "mouse_T_gamma"
            if (recomb_type is None):
                recomb_type = None # ??? Not sure what is teh most appropriate model
            cdr3s = self.clone_df.cdr3_g_aa

            if not cdr3_only:
                v_genes = self.clone_df.v_g_gene
                j_genes = self.clone_df.j_g_gene
            else:
                v_genes = None
                j_genes = None

        if chain  ==  "delta":
            raise ValueError("SORRY:OLGA default files do not yet support delta TCRs")
            if (chain_folder is None):
                if (chain_folder is None):
                    if self.organism is 'human':
                        chain_folder = "human_T_delta"
                    elif self.organism is 'mouse':
                        chain_folder = "mouse_T_delta"
            if (recomb_type is None):
                recomb_type = None # ??? Not sure what is teh most appropriate model
            cdr3s = self.clone_df.cdr3_d_aa

            if not cdr3_only:
                v_genes = self.clone_df.v_d_gene
                j_genes = self.clone_df.j_d_gene
            else:
                v_genes = None
                j_genes = None


        # initializes the appropriate olga genomic model
        my_olga_model = pgen.OlgaModel(chain_folder = chain_folder,
                                       recomb_type = recomb_type)
        # computes pgen from clone_df
        olga_pgens = my_olga_model.compute_aa_cdr3_pgens(cdr3s,
                                                         v_genes,
                                                         j_genes)

        if chain is "alpha":
            self.clone_df['cdr3_a_aa_pgen'] = pd.Series(olga_pgens)
        if chain is "beta":
            self.clone_df['cdr3_b_aa_pgen'] = pd.Series(olga_pgens)
        if chain is "gamma":
            self.clone_df['cdr3_g_aa_pgen'] = pd.Series(olga_pgens)
        if chain is "delta":
            self.clone_df['cdr3_d_aa_pgen'] = pd.Series(olga_pgens)

        return(pd.Series(olga_pgens))
    
    def archive(self, dest = "default_archive", dest_tar_name = "default_archive.tar.gz" ):
        """ Use Zipdist2 to Make an Archive.tar.gz """
        self.cell_df_index = self.cell_df.index.copy()
        self.cell_df = self.cell_df.reset_index()
        z = Zipdist2(name = dest_tar_name , target = self)
        z._save(dest = dest, dest_tar = dest_tar_name )
        sys.stdout.write(f"\tArchiving your TCRrep using Zipdist2 in [{dest_tar_name}]\n")
    
    def rebuild(self,  dest_tar_name = "default_archive.tar.gz" ):
        """ Use Zipdist2 to Make an Archive.tar.gz """
        #tr = TCRrep(cell_df=df.iloc[0:0,:], chains=chains, organism='mouse')
        z = Zipdist2(name = "default_archive", target = self)
        z._build(dest_tar = dest_tar_name , target = self)
        
        # VALIDATION OF INPUTS
        # check that chains are valid.
        self._validate_organism()
        self._validate_chains()
        # check that  is a pd.DataFrame
        self._validate_cell_df()

        # RE INIT the REFERENCE DB see repertoire_db.py
        self.generate_ref_genes_from_db(self.db_file)
       


    def tcrdist2(self, metric = "nw", processes = None, weights = None, dump = False, reduce = True, save = False, dest = "default_archive", dest_tar_name = "default_archive.tar.gz" ):
        """
        Automated calculation of single chain and paired chain tcr-distances

        Parameters
        ----------
        metric : str
            specified metric, currently only "nw" and "hamming" are supported 
            (see notes for legacy methods)
        processes : int
            number of cpus to use; the default is greedy and will use half of available
        weights : dict
            override cdr weightings
        dump : bool
            if True, dump intermediate cdr1, cdr2, and pmhc pairwise matrices
        reduce : bool
            if True, converts distance matrices to a smaller data type.
        save : bool
            if True, saves intermediate files to dest
        dest : str
            path to save components
        

        
        Notes
        -----
        tcrdist2 is a method to help new-users run tcrdist2 with sensible defaults.

        Distance metrics are highly customizable. 
        Consult the `docs <https://tcrdist2.readthedocs.io>`_ for more information.
        
        To compute Dash et al. 2017 style tcrdistance, instead of tcrdist2,
        use commands:
        TCRrep._tcrdist_legacy_method_alpha_beta()
        TCRrep._tcrdist_legacy_method_beta()
        TCRrep._tcrdist_legacy_method_alpha()
        TCRrep._tcrdist_legacy_method_gamma_delta()
        TCRrep._tcrdist_legacy_method_gamma()
        TCRrep._tcrdist_legacy_method_delta()
        """

        # Default to use all available processes
        if processes is None:
            max_threads = multiprocessing.cpu_count() 
            processes = max_threads // 2
            sys.stdout.write(f"trcdist2 detected {max_threads } available cpus/threads.\n")
            sys.stdout.write(f"\tTCRrep use parallel processing, setting default to use {processes} cpus/threads.\n")
            sys.stdout.write(f"\tThe `processes` arg of TCRrep.tcrdist2() can be set manually\n")

        for chain in self.chains:
            self.infer_cdrs_from_v_gene(chain=chain,  imgt_aligned=True)
        
        if weights is None:
            weights = {'cdr1_a_aa':1,
                        'cdr2_a_aa':1,
                        'cdr3_a_aa':3,
                        'pmhc_a_aa':1,
                        'cdr1_b_aa':1,
                        'cdr2_b_aa':1,
                        'cdr3_b_aa':3,
                        'pmhc_b_aa':1,
                        'cdr1_g_aa':1,
                        'cdr2_g_aa':1,
                        'cdr3_g_aa':3,
                        'pmhc_g_aa':1,
                        'cdr1_d_aa':1,
                        'cdr2_d_aa':1,
                        'cdr3_d_aa':3,
                        'pmhc_d_aa':1,
                        'v_a_gene':0, 
                        'j_a_gene':0,
                        'v_b_gene':0, 
                        'j_b_gene':0,
                        'v_g_gene':0, 
                        'j_g_gene':0,
                        'v_d_gene':0, 
                        'j_d_gene':0,
                        'cdr3_a_nucseq':0,
                        'cdr3_b_nucseq':0,
                        'cdr3_g_nucseq':0,
                        'cdr3_d_nucseq':0}

        index_cdrs = [k for k in weights.keys() if k in self.cell_df.columns]

        for x in ['clone_id', 'subject', 'epitope']:
            assert 'clone_id' in self.cell_df.columns, f"{x} must be in TCRrep.cell_df"

        self.index_cols = ['clone_id', 'subject', 'epitope'] + index_cdrs
        sys.stdout.write("Deduplicating your TCRrep.cell_df to make TCRrep.clone_df.\n")

        self.deduplicate()
        sys.stdout.write(f"Computing pairwise matrices for multiple Complementarity Determining Regions (CDRs):.\n")
        for chain in self.chains:
            sys.stdout.write(f"\tComputing pairwise matrices for cdrs within the {chain}-chain using the {metric} metric.\n")
            self.compute_pairwise_all(chain = chain, metric = metric, processes = processes)
        sys.stdout.write("Calculating composite tcrdistance measures:\n")
        self.compute_paired_tcrdist( chains=self.chains, store_result=False)
        for chain in self.chains:
            sys.stdout.write(f"\tSingle chain pairwise tcrdistances are in attribute : TCRrep.pw_{chain}\n")
        sys.stdout.write(f"\tCombined pairwise tcrdistances are in attribute     : TCRrep.pw_tcrdist\n")
        sys.stdout.write(f"\tCDR specific tcrdistances are in attributes, e.g.,  : TCRrep.cdr3_{chain[0]}_aa_pw\n")


        # <dump> boolean controls whether we dump easy to recalculate cdr1, cdr2, pmhc 

        # <shrink> boolean controls whether we convert distance matrices 
        # to a smaller data type.
        if reduce:
            data_type = 'int16'
            sys.stdout.write(f"Reducing File Size: `reduce` argumment set to {reduce}:\n")
            self.reduce_file_size( data_type = data_type)
              
        # pairwise matices, which most users will never again.
        if dump:
            sys.stdout.write(f"Cleanup: `dump` argument set to {dump}. Dumping individual CDR specific distance matrices:\n")
            for i in index_cdrs:
                if i.startswith("cdr1") or i.startswith("cdr2") or i.startswith("pmhc"):
                    if i.endswith("aa"):
                        i = f"{i}_pw"
                        sys.stdout.write(f"\tDumping : {i}\n")                    
                        self.__dict__[i] = None
        if save:
            sys.stdout.write(f"Archiving your TCRrep using Zipdist2 (save = {save})\n")
            # To avoid = ValueError: feather does not support serializing a non-default index for the index; you can .reset_index() to make the index into column(s)
            self.archive(dest = dest, dest_tar_name = dest_tar_name )
            sys.stdout.write(f"\tArchiving your TCRrep using Zipdist2 in [{dest_tar_name}]\n")


        sys.stdout.write(f"TCRrep.tcrdist2() COMPLETED SUCCESSFULLY, see the docs for Analysis steps!\n")

        

            




    def compute_pairwise_all(self,
                             chain,
                             compute_specific_region = None,
                             metric = "hamming",
                             processes = 2,
                             user_function = None,
                             to_matrix = True,
                             **kwargs):
        """
        Computes pairwise distances for all regions on a given
        chain or for a specific region on that chain.

        Parameters
    	----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        compute_specific_region : string
            optional string (e.g. "cdr2_a_aa") to over-ride function behavior
            and compute only a single region
        metric : string
            'nw', 'hamming', or 'custom' (or if legacy tcrdist is to be calculated,
            "tcrdist_cdr3", "tcrdist_cdr1", "tcrdist_cdr2",
            "tcrdist_cdr2.5", "tcrdist_pmhc" can be supplied. WARNING:
            imgt_aligned must be set to True in tr.infer_cdrs_from_v_gene().
        processes : int
            int for number of available cpu for multiprocessing (to see available
            try multiprocessing.cpu_count())
        user_function : function
            function for a custom distance metric on two strings (This is
            an advanced option, so don't use this unless you are absolutely
            sure what you are doing; metric arg must be set to 'custom').
        to_matrix : boolean
            True will return pairwise distance as result as a 2D ndarray




        Notes
    	-----

        Uses _assign_pw_result to assign self.[cdr3|cdr2|cdr1|pmhc]_[a|b|d|g]_aa_pw objects


        Examples
    	--------
        >>> testrep = TCRrep(cell_df = example_df, organism = "human", chains= ["alpha","beta"])
        >>> testrep.infer_cdrs_from_v_gene(chain = "alpha")
        >>> testrep.infer_cdrs_from_v_gene(chain = "beta")
        >>> testrep.index_cols = testrep.index_cols + ['cdr1_a_aa','cdr2_a_aa','pmhc_a_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']
        >>> testrep.deduplicate()
        >>> testrep.compute_pairwise_all(chain = "alpha", metric= "hamming")
        >>> testrep.compute_pairwise_all(chain = "beta", metric= "hamming")

        alternatively, compute each region one by one

        >>> testrep.compute_pairwise_all(chain = "beta", compute_specific_region="cdr1_b_aa")
        >>> testrep.compute_pairwise_all(chain = "alpha", compute_specific_region="cdr2_a_aa")

        """

        # validate chain argument passed
        self._validate_chain(chain)

        if metric in ["tcrdist_cdr3", "tcrdist_cdr1", "tcrdist_cdr2",
                      "tcrdist_cdr2.5", "tcrdist_pmhc"]:
            if not self.imgt_aligned_status:
                raise ValueError("imgt_aligned must be set to True in tr.infer_cdrs_from_v_gene()")

        # If compute_specific_region is None, then the behavior is to loop through the a list regions.
        if compute_specific_region is None:
            index_col_from_chain = {'alpha' : ['cdr3_a_aa', 'cdr2_a_aa',
                                               'cdr1_a_aa', 'pmhc_a_aa'],
                                    'beta'  : ['cdr3_b_aa', 'cdr2_b_aa',
                                               'cdr1_b_aa', 'pmhc_b_aa'],
                                    'gamma' : ['cdr3_g_aa', 'cdr2_g_aa',
                                               'cdr1_g_aa', 'pmhc_g_aa'],
                                    'delta' : ['cdr3_d_aa', 'cdr2_d_aa',
                                               'cdr1_d_aa', 'pmhc_d_aa']}
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
                               replacement_weights = {},
                               store_result = False):
        """
        Computes tcrdistance metric combining distances metrics across multiple
        T Cell Receptor CDR regions.

        Parameters
        ----------
        chains : list
            list of strings containing some combination of 'alpha', 'beta',
            'gamma', and 'delta'
        replacement_weights : dictionary
            optional dictionary of the form {'cdr1_a_aa_pw':1, 'cdr2_a_aa_pw':1}
            used to place greater weight on certain TCR regions. The default
            is a weight of 1.
        store_result : boolean
            True will store results to
            :py:attr:`TCRrep.stored_tcrdist`

        Returns
        -------
        r : dictionary
            a dictionary with keys paired_tcrdist points to a 2D
            tcrdist np.ndarray and paired_tcrdist_weights pointing to
            dictionary of weights. See notes.

        Notes
        -----

        Calling this function assigns results to
        `TCRrep.paired_tcrdist` and
        `TCRrep.paired_tcrdist_weights`
        and stores r to
        `TCRrep.stored_tcrdist`

        In addition it returns a dictionary with keys `paired_tcrdist` 2D
        tcrdist np.array and `paired_tcrdist_weights`
        a dictionary of regions and relative weights:

        {'paired_tcrdist': array([[ 0., 76., 80.,..., 89., 89., 87.],
                                [ 76., 0., 60., ..., 81., 75., 43.],
                                [ 80., 60., 0., ..., 59., 81., 77.],
                                ...,
                                [ 89., 81., 59.,  ..., 0., 60., 58.],
                                [ 89., 75., 81.,   ..., 60., 0., 40.],
                                [ 87., 43., 77., ..., 58., 40., 0.]]),
        'paired_tcrdist_weights': {'cdr1_a_aa_pw': 1,
                                   'cdr1_b_aa_pw': 2,
                                   'cdr2_a_aa_pw': 1,
                                   'cdr2_b_aa_pw': 2,
                                   'cdr3_a_aa_pw': 2,
                                   'cdr3_b_aa_pw': 4,
                                   'pmhc_a_aa_pw': 1,
                                   'pmhc_b_aa_pw': 2}}

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

        alpha_keys = [k for k in list(weights.keys()) if k.endswith("a_aa_pw")]
        beta_keys  = [k for k in list(weights.keys()) if k.endswith("b_aa_pw")]
        gamma_keys = [k for k in list(weights.keys()) if k.endswith("g_aa_pw")]
        delta_keys = [k for k in list(weights.keys()) if k.endswith("d_aa_pw")]

        # for single chain computation, results in TCRrep.pw_alpha, TCRrep.pw_beta, TCRrep.pw_gamma, and or TCRrep.pw_delta, 
        if 'alpha' in chains:
            tcrdist = np.zeros(self.cdr3_a_aa_pw.shape)
            for k in alpha_keys:
                try:
                    tcrdist = self.__dict__[k]*weights[k] + tcrdist
                except KeyError:
                    warnings.warn("tcrdist was calculated without: '{}' because pairwise distances haven't been computed for this region:".format(k))
            self.pw_alpha = tcrdist
        if 'beta' in chains:
            tcrdist = np.zeros(self.cdr3_b_aa_pw.shape)
            for k in beta_keys:
                try:
                    tcrdist = self.__dict__[k]*weights[k] + tcrdist
                except KeyError:
                    warnings.warn("tcrdist was calculated without: '{}' because pairwise distances haven't been computed for this region:".format(k))
            self.pw_beta = tcrdist
        if 'gamma' in chains:
            tcrdist = np.zeros(self.cdr3_g_aa_pw.shape)
            for k in gamma_keys:
                try:
                    tcrdist = self.__dict__[k]*weights[k] + tcrdist
                except KeyError:
                    warnings.warn("tcrdist was calculated without: '{}' because pairwise distances haven't been computed for this region:".format(k))
            self.pw_gamma = tcrdist
        if 'delta' in chains:
            tcrdist = np.zeros(self.cdr3_d_aa_pw.shape)
            for k in delta_keys:
                try:
                    tcrdist = self.__dict__[k]*weights[k] + tcrdist
                except KeyError:
                    warnings.warn("tcrdist was calculated without: '{}' because pairwise distances haven't been computed for this region:".format(k))
            self.pw_delta = tcrdist
        
        # For combined chain tcrdist, restults in TCRrep.paired_tcrdist and TCRrep.pw_tcrdist
        full_keys = []
        if 'alpha' in chains:
            full_keys = full_keys + alpha_keys
        if 'beta' in chains:
            full_keys = full_keys + beta_keys
        if 'gamma' in chains:
            full_keys = full_keys + gamma_keys
        if 'delta' in chains:
            full_keys = full_keys + delta_keys

        # initialize tcrdist matrix size
        for k in full_keys:
            try:
                tcrdist = np.zeros(self.__dict__[k].shape)
                break
            except KeyError:
                pass

        for k in full_keys:
            try:
                tcrdist = self.__dict__[k]*weights[k] + tcrdist
            except KeyError:
                warnings.warn("tcrdist was calculated without: '{}' because pairwise distances haven't been computed for this region:".format(k))
                pass
     
        # keep 'paired_tcrdist' to avoid breaking tests
        self.paired_tcrdist = tcrdist
        self.pw_tcrdist = tcrdist
        self.paired_tcrdist_weights = {k:weights[k] for k in full_keys}
        
        # Typically we don't want to store different tcrdistance in the same repertoire, but
        r = {'paired_tcrdist' : tcrdist,
            'paired_tcrdist_weights' : {k:weights[k] for k in full_keys}}
        if store_result:
            self.stored_tcrdist.append(r)
        return(r)

    def compute_pairwise(self,
                         chain,
                         metric = "nw",
                         processes = 2,
                         user_function = None,
                         to_matrix = True,
                         **kwargs):
        """
        Early Function to be replaced with compute_pairwise_all.
        TODO: Rewrite test and remove.
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
    
    def generate_cluster_index(self, t = 75, criterion = "distance", method =  "complete", append_counts = False):
        """
        Add 'cluster_index' column to TCRrep.clone_df 

        Parameters
        ----------

        t : int
            scipy.cluster.hierarchy.fcluster param t
        criterion : str 
            scipy.cluster.hierarchy.fcluster param criterion 
        method : str 
            scipy.cluster.linkage parma method 
        
        Notes
        -----
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html
        """
        compressed_dmat = scipy.spatial.distance.squareform(self.paired_tcrdist, force = "vector")
        Z = linkage(compressed_dmat, method = "complete")
        cluster_index = fcluster(Z, t = t, criterion = criterion)
        assert len(cluster_index) == self.clone_df.shape[0]
        assert len(cluster_index) == self.paired_tcrdist.shape[0]
        self.clone_df['cluster_index'] = cluster_index
        if append_counts:
            self._append_cluster_count()
            self._append_seq_counts_per_cluster()

    def _append_cluster_count(self):
        """
        Appends the number of clones in a cluster to each row of TCRrep.clone_df
        """
        cluster_count = self.clone_df.cluster_index.value_counts().\
            reset_index().\
            rename(columns = {'index':'cluster_index', "cluster_index": "cluster_count"}).\
            copy()
        self.clone_df = self.clone_df.merge(cluster_count, how= "left", left_on = "cluster_index", right_on = "cluster_index")

    def _append_seq_counts_per_cluster(self):
        """
        Appends the sum of seq counts per cluster to each row of TCRrep.clone_df
        """
        seq_counts = self.clone_df.\
            groupby(['cluster_index'])['count'].sum().\
            reset_index().\
            rename(columns = {'count':'seq_count'})

        self.clone_df = self.clone_df.merge(seq_counts, how = "left", left_on = "cluster_index", right_on = "cluster_index")    

    def ispublic(self, gr, var = "subject", n = 1):
        """
        Return True if a cluster public, defined as comprised of members from multiple individuals 
        or cell subsets (e.g., CD4/CD8) 
        
        Parameters
        ----------
        gr : group
            within pandas Data.Frame.groupby
        var : str
            variable name of class that group most transcend to be considered public
        m : int
            number of unique values of selected variable to be considered public
        
        Returns
        -------
        r : bool
            True if a cluster public
        
        """
        r = len(gr[var].value_counts()) > n
        if r:
            return 'public'
        else:
            return 'private'

    def get_func_stat(self, gr, var, func = np.median, **kwargs):
        """
        get summary statistic by applying a func to a group
        
        Parameter
        ---------
        gr : group
            within pandas Data.Frame.groupby
        var : str
            variable name of class that group most transcend to be considered public
        func : function
            function that can operate on a series or list of values specified by var
            
        Returns
        -------
        r : float or int
        
        """
        r = func(gr[var], **kwargs)
        return r

    def get_cluster_summary(self, df=None, groupvar = 'cluster_index'):
        """get_cluster_pgen_and_count_summary """
        if df is None:
            df = self.clone_df.copy()
        cluster_summary = list()
        assert groupvar in df.columns
        assert "pgen" in df.columns
        assert "count" in df.columns
        for name, group in df.groupby([groupvar]):
            public = self.ispublic(group, "subject")
            cluster_summary.append({"cluster_index" : name, 
                                "public"        : public,
                                "min_pgen"      : self.get_func_stat(gr = group, var = "pgen", func = np.min),
                                "median_pgen"   : self.get_func_stat(gr = group, var = "pgen", func = np.median),
                                "max_pgen"      : self.get_func_stat(gr = group, var = "pgen", func = np.max), 
                                "cluster_count"  : group.shape[0],
                                "seq_count"     : self.get_func_stat(gr = group, var = "count", func = np.sum),
                                "seq_min"       : self.get_func_stat(gr = group, var = "count", func = np.min),
                                "seq_median"    : self.get_func_stat(gr = group, var = "count", func = np.median),
                                "seq_max"       : self.get_func_stat(gr = group, var = "count", func = np.max)})

        cluster_summary = pd.DataFrame(cluster_summary)
        self.cluster_summary = cluster_summary
        return cluster_summary

    def tsne(self, X = None, n_components=2 , random_state = 310, axis_names = ["tSNE1","tSNE2"]):
        warnings.warn("RUNNING sklearn.manifold.TSNE WHICH MAY TAKE A FEW MINUTES")
        from sklearn.manifold import TSNE
        if X is None:
            X = self.paired_tcrdist        
<<<<<<< HEAD
        X_embedded = TSNE(n_components=n_components, metric = 'precomputed', random_state = random_state).fit_transform(X)
=======
        X_embedded = TSNE(n_components=n_components, metric ='precomputed', random_state = random_state).fit_transform(X)
>>>>>>> e6a232884c70488d83e6b543169eb6451a9becd9
        tsne_df = pd.DataFrame(X_embedded, columns = axis_names )
        assert(tsne_df.shape[0] == self.clone_df.shape[0])
        self.clone_df = pd.concat([self.clone_df, tsne_df], axis = 1)
    
    def mds(self, X = None, n_components=2 , dissimilarity='precomputed', axis_names = ["MDS1","MDS2"]):
        warnings.warn("RUNNING sklearn.manifold.MDS WHICH MAY TAKE A FEW MINUTES")
        from sklearn.manifold import MDS
        if X is None:
            X = self.paired_tcrdist     
        X_embedded_mds = MDS(n_components=n_components, dissimilarity=dissimilarity).fit_transform(X)
        mds_df = pd.DataFrame(X_embedded_mds, columns = axis_names)
        assert(mds_df.shape[0] == self.clone_df.shape[0])
        self.clone_df = pd.concat([self.clone_df, mds_df], axis = 1)

    def _validate_organism(self):
        if self.organism not in ["human", "mouse"]:
            raise ValueError("organism must be 'mouse' or 'human'")

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
            self.cdr3_a_aa_smat = 'blosum62'
            self.cdr2_a_aa_smat = 'blosum62'
            self.cdr1_a_aa_smat = 'blosum62'
            self.pmhc_a_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_a_aa")

        if 'beta' in self.chains:
            self.cdr3_b_aa_smat = 'blosum62'
            self.cdr2_b_aa_smat = 'blosum62'
            self.cdr1_b_aa_smat = 'blosum62'
            self.pmhc_b_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_b_aa")

        if 'gamma' in self.chains:
            self.cdr3_g_aa_smat = 'blosum62'
            self.cdr2_g_aa_smat = 'blosum62'
            self.cdr1_g_aa_smat = 'blosum62'
            self.pmhc_g_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_g_aa")

        if 'delta' in self.chains:
            self.cdr3_d_aa_smat = 'blosum62'
            self.cdr2_d_aa_smat = 'blosum62'
            self.cdr1_d_aa_smat = 'blosum62'
            self.pmhc_d_aa_smat = 'blosum62'
            self.index_cols.append("cdr3_d_aa")



    def _get_smat(self, chain, index_col):
        """
        Gets the correct substitution matrix (smat) based on chain and column

        Parameters
        ----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        index_col : string
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
                smat = 'blosum62'
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
                smat = 'blosum62'
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
                smat = 'blosum62'
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
                smat = 'blosum62'
                warnings.warn("Using default parasail.blosum62 because chain: '{}' does not matches region: '{}'".format(index_col, chain, index_col))

        return(smat)


    def _assign_pw_result(self, pw, chain, index_col):
        """
        Assigns pairwise result to TCRrep attribute based on chain and index_col

        Parameters
        ----------
        chain : string
            'alpha', 'beta', 'gamma', or 'delta'
        index_col : string
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

    def _drop_smats(self):
        """
        Need to drop ctypes if you are to pickle or copy this instance
        """
        smats = [ k for k in self.__dir__() if k.endswith("aa_smat") ]
        for k in smats:
            self.__dict__[k] = None

    def _pickle(self, filename):
        self._drop_smats()
        pickle.dump(self,  open(filename , "wb") )
        warnings.warn("all smats dropped because they are C objects that can't be pickled. reassign with _initialize_chain_specific_attributes()")

    def _tcrdist_legacy_method_alpha_beta(self, processes = 1):
        """
        Runs the legacy tcrdist pairwise comparison

        Arguments
        ---------
        processes : int


        Notes
        -----
        # CALCULATE tcrdist distance metric. Here we show all the manual steps to
        # implement the original Dash et al. tcrdist approach.

        # To do this we calculate distance for each CDR separately, and
        # we use the metric "tcrdist_cdr3" for the cdr3 and "tcrdist_cdr1"
        # everywhere else
        """


        if "gamma" in self.chains or "delta" in self.chains:
            raise ValueError("TCRrep.chains contains `gamma`. You might want "\
            "TCRrep._tcrdist_legacy_method_gamma_delta")

        self.compute_pairwise_all(chain = "alpha",                        # <11
                                 metric = 'tcrdist_cdr3',
                                 compute_specific_region = 'cdr3_a_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "alpha",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_a_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "alpha",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_a_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "alpha",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_a_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = 'tcrdist_cdr3',
                                 #user_function = tcrdist_metric_align_cdr3s_false,
                                 compute_specific_region = 'cdr3_b_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_b_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_b_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_b_aa',
                                 processes = processes)

        distA = self.compute_paired_tcrdist(replacement_weights= {'cdr3_a_aa_pw': 1,
                                                                'cdr2_a_aa_pw': 1,
                                                                'cdr1_a_aa_pw': 1,
                                                                'pmhc_a_aa_pw': 1,
                                                                'cdr3_b_aa_pw': 0,
                                                                'cdr2_b_aa_pw': 0,
                                                                'cdr1_b_aa_pw': 0,
                                                                'pmhc_b_aa_pw': 0},
                                                                 chains = ["alpha", "beta"])['paired_tcrdist'].copy()

        distB = self.compute_paired_tcrdist(replacement_weights= {'cdr3_a_aa_pw': 0,
                                                                'cdr2_a_aa_pw': 0,
                                                                'cdr1_a_aa_pw': 0,
                                                                'pmhc_a_aa_pw': 0,
                                                                'cdr3_b_aa_pw': 1,
                                                                'cdr2_b_aa_pw': 1,
                                                                'cdr1_b_aa_pw': 1,
                                                                'pmhc_b_aa_pw': 1},
                                                                chains = ["alpha", "beta"])['paired_tcrdist'].copy()

        # Calling tr.compute_paired_tcrdist() computs the
        # the final paired chain TCR-distance which is stored as
        # tr.paired_tcrdist, which we confirm is simply the sum of distA and distB
        self.compute_paired_tcrdist()
        assert np.all(((distA + distB) - self.paired_tcrdist) == 0)

        # tr.paired_tcrdist and distA, distB are np arrays, but we will want to work with as a pandas DataFrames
        self.dist_a = pd.DataFrame(distA, index = self.clone_df.clone_id, columns = self.clone_df.clone_id)
        self.dist_b = pd.DataFrame(distB, index = self.clone_df.clone_id, columns = self.clone_df.clone_id)

    def _tcrdist_legacy_method_alpha(self, processes = 1):
        """
        Runs the legacy tcrdist pairwise comparison

        Arguments
        ---------
        processes : int


        Notes
        -----
        # CALCULATE tcrdist distance metric. Here we show all the manual steps to
        # implement the original Dash et al. tcrdist approach.

        # To do this we calculate distance for each CDR separately, and
        # we use the metric "tcrdist_cdr3" for the cdr3 and "tcrdist_cdr1"
        # everywhere else
        """


        if "gamma" in self.chains or "delta" in self.chains:
            raise ValueError("TCRrep.chains contains `gamma`. You might want "\
            "TCRrep._tcrdist_legacy_method_gamma_delta")

        self.compute_pairwise_all(chain = "alpha",                        # <11
                                    metric = 'tcrdist_cdr3',
                                    compute_specific_region = 'cdr3_a_aa',
                                    processes = processes)
        self.compute_pairwise_all(chain = "alpha",                        # 11
                                    metric = "tcrdist_cdr1",
                                    compute_specific_region = 'cdr1_a_aa',
                                    processes = processes)
        self.compute_pairwise_all(chain = "alpha",                        # 11
                                    metric = "tcrdist_cdr1",
                                    compute_specific_region = 'cdr2_a_aa',
                                    processes = processes)
        self.compute_pairwise_all(chain = "alpha",                        # 11
                                    metric = "tcrdist_cdr1",
                                    compute_specific_region = 'pmhc_a_aa',
                                    processes = processes)

        distA = self.compute_paired_tcrdist(replacement_weights= {'cdr3_a_aa_pw': 1,
                                                                'cdr2_a_aa_pw': 1,
                                                                'cdr1_a_aa_pw': 1,
                                                                'pmhc_a_aa_pw': 1},
                                                                chains = ["alpha"])['paired_tcrdist'].copy()
                                                                        # Calling tr.compute_paired_tcrdist() computs the
        # the final paired chain TCR-distance which is stored as
        # tr.paired_tcrdist, which we confirm is simply the sum of distA and distB
        self.compute_paired_tcrdist()
        assert np.all((distA - self.paired_tcrdist) == 0)

        # tr.paired_tcrdist and distA, distB are np arrays, but we will want to work with as a pandas DataFrames
        self.dist_a = pd.DataFrame(distA, index = self.clone_df.clone_id, columns = self.clone_df.clone_id)


    def _tcrdist_legacy_method_beta(self, processes = 1):
        """
        Runs the legacy tcrdist pairwise comparison

        Arguments
        ---------
        processes : int


        Notes
        -----
        # CALCULATE tcrdist distance metric. Here we show all the manual steps to
        # implement the original Dash et al. tcrdist approach.

        # To do this we calculate distance for each CDR separately, and
        # we use the metric "tcrdist_cdr3" for the cdr3 and "tcrdist_cdr1"
        # everywhere else
        """


        if "gamma" in self.chains or "delta" in self.chains:
            raise ValueError("TCRrep.chains contains `gamma`. You might want "\
            "TCRrep._tcrdist_legacy_method_gamma_delta")

        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = 'tcrdist_cdr3',
                                 #user_function = tcrdist_metric_align_cdr3s_false,
                                 compute_specific_region = 'cdr3_b_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_b_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_b_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "beta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_b_aa',
                                 processes = processes)

        distB = self.compute_paired_tcrdist(replacement_weights= {
                                                                'cdr3_b_aa_pw': 1,
                                                                'cdr2_b_aa_pw': 1,
                                                                'cdr1_b_aa_pw': 1,
                                                                'pmhc_b_aa_pw': 1},
                                                                chains = ["beta"])['paired_tcrdist'].copy()

        self.compute_paired_tcrdist()
        assert np.all(((distB) - self.paired_tcrdist) == 0)

        # tr.paired_tcrdist and distA, distB are np arrays, but we will want to work with as a pandas DataFrames
        self.dist_b = pd.DataFrame(distB, index = self.clone_df.clone_id, columns = self.clone_df.clone_id)



    def _tcrdist_legacy_method_gamma_delta(self, processes = 1):
        """
        Runs the legacy tcrdist pairwise comparison gamma/delta

        Arguments
        ---------
        processes : int


        Notes
        -----
        # CALCULATE tcrdist distance metric. Here we show all the manual steps to
        # implement the original Dash et al. tcrdist approach.

        # To do this we calculate distance for each CDR separately, and
        # we use the metric "tcrdist_cdr3" for the cdr3 and "tcrdist_cdr1"
        # everywhere else
        """
        self.compute_pairwise_all(chain = "gamma",                        # <11
                                 metric = 'tcrdist_cdr3',
                                 compute_specific_region = 'cdr3_g_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "gamma",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_g_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "gamma",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_g_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "gamma",                        # 11
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_g_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "delta",                         # 12
                                 metric = 'tcrdist_cdr3',
                                 compute_specific_region = 'cdr3_d_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "delta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr1_d_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "delta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'cdr2_d_aa',
                                 processes = processes)
        self.compute_pairwise_all(chain = "delta",                         # 12
                                 metric = "tcrdist_cdr1",
                                 compute_specific_region = 'pmhc_d_aa',
                                 processes = processes)

        distA = self.compute_paired_tcrdist(replacement_weights= {'cdr3_g_aa_pw': 1,
                                                                'cdr2_g_aa_pw': 1,
                                                                'cdr1_g_aa_pw': 1,
                                                                'pmhc_g_aa_pw': 1,
                                                                'cdr3_d_aa_pw': 0,
                                                                'cdr2_d_aa_pw': 0,
                                                                'cdr1_d_aa_pw': 0,
                                                                'pmhc_d_aa_pw': 0},
                                                                chains = ["gamma", "delta"])['paired_tcrdist'].copy()

        distB = self.compute_paired_tcrdist(replacement_weights= {'cdr3_g_aa_pw': 0,
                                                                'cdr2_g_aa_pw': 0,
                                                                'cdr1_g_aa_pw': 0,
                                                                'pmhc_g_aa_pw': 0,
                                                                'cdr3_d_aa_pw': 1,
                                                                'cdr2_d_aa_pw': 1,
                                                                'cdr1_d_aa_pw': 1,
                                                                'pmhc_d_aa_pw': 1},
                                                                chains = ["gamma", "delta"] )['paired_tcrdist'].copy()

        # Calling tr.compute_paired_tcrdist() computs the
        # the final paired chain TCR-distance which is stored as
        # tr.paired_tcrdist, which we confirm is simply the sum of distA and distB
        self.compute_paired_tcrdist(chains = ["gamma", "delta"])
        assert np.all(((distA + distB) - self.paired_tcrdist) == 0)

        # tr.paired_tcrdist and distA, distB are np arrays, but we will want to work with as a pandas DataFrames
        self.dist_g = pd.DataFrame(distA, index = self.clone_df.clone_id, columns = self.clone_df.clone_id)
        self.dist_d = pd.DataFrame(distB, index = self.clone_df.clone_id, columns = self.clone_df.clone_id)    

    def drop_stored_results(self):
        """
        Sets list of stored results to None.
        This will not effect the individual CDR pair-wise matrices or
        the most recent self.paired_tcrdist

        """
        self.stored_tcrdist = None

    def reduce_file_size(self, attributes = None, data_type = 'int16'):
        """
        Parameters
        ----------

        attributes : list
            The list of numpy attributes to change data type.
            If left blank all numpy attributes will be changes
        data_type : str
            numpy data type default = 'int16' which will accommodate
            tcrdistances but save considerable space in memory.

        Notes
        -----
        This is useful for reducing the file size and memory usage of the
        large matrices stored as attributes of the TCRrep class. E.g.,
        a 10,000 x 10,000 pairwise matrix can be shrunk from > 1GB to
        300 MB by converting float to int16.

        Numpy supports a much greater variety of numerical types than Python
        does. Some are:

        int16	Integer (-32768 to 32767)
        int32	Integer (-2147483648 to 2147483647)

        See https://docs.scipy.org/doc/numpy-1.10.0/user/basics.types.html for
        details on other types.

        """
        # If attributes are not user supplied, all TCRrep attributes will be
        # considered for type conversion.

        if attributes is None:
            attributes = list(self.__dict__.keys())

        if not isinstance(attributes, list):
            raise TypeError("< attributes > argument must be supplied as a list")

        # However, only np.arrays will be transformed to the new data type
        for k in attributes:
            if isinstance(getattr(self,k), np.ndarray):
                sys.stdout.write("\tReducing : {} to {}.\n".format(k, data_type))
                self._reduce_file_size_of_np_attribute(k, data_type)


    def _reduce_file_size_of_np_attribute(self, attribute,  data_type = 'int16'):
        """
        Parameters
        ----------
        attribute : str
            TCRrep attribute
        data_type : str
            numpy datatype
        """
        x = getattr(self, attribute).astype(data_type)
        setattr(self, attribute, x)

    def save_as_hdf5(self, hdf5_file):
        """
        Save class attributes as a hdf5 file

        Parameters
        ----------
        hdf5_file : string
            specifies the hdf5 file to used to save TCRrep attribute
        """
        hdf = pd.HDFStore(hdf5_file)
        all_attrs = self.__dict__.keys()
        all_types = [type(getattr(self,x)) for x in all_attrs]

        attributes_dict = {k : getattr(self,k) for k,t in zip(all_attrs, all_types) if t in [bool, str, int, list, dict]}
        attributes_dict['attributes'] = list(all_attrs)

        for attr, ty in zip(all_attrs, all_types):
            x = getattr(self, attr)
            if not isinstance(x, parasail.bindings_v2.Matrix):
                if isinstance(x, pd.core.frame.DataFrame):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        hdf.put(attr, x, format='table', data_columns=True)
                if isinstance(x, np.ndarray):
                    if 'pw' in attr:
                        out = scipy.spatial.distance.squareform(x, force='tovector', checks=False)
                        out = pd.DataFrame(out)
                    else:
                        out = pd.DataFrame(x)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        hdf.put(attr, out, format='table')
        hdf.get_storer('cell_df').attrs.metadata = json.dumps(attributes_dict)
        hdf.close()

    def rebuild_from_hdf5(self, hdf5_file, verbose=False):
        """
        Use hdf5 file to rebuild a TCRrep instance.

        Parameters
        ----------
        hdf5_file : string
            specifies the hdf5 file used to repopulate a TCRdist object

        Notes
        -----
        pandas and numpy objects are stored as datasets within the hdf5 file.

        str, list, bool, and float type attributes are stored in a json
        file that is set as metadata for the archived cell_df dataset.
        Within this metadata, the attributes contains a list of all class
        attributes.
        """
        # Currently these are the only acceptable attributes that can be
        # read from the HDF5
        attr_to_type = {
        'db_file' : str,
        'cell_df': pd.core.frame.DataFrame,
        'chains': list,
        'organism': str,
        'pwdist_df': None,
        'clone_df': pd.core.frame.DataFrame,
        'index_cols': list,
        'stored_tcrdist': None,
        'paired_tcrdist': np.ndarray,
        'paired_tcrdist_weights': dict,
        'meta_cols': None,
        'project_id': str,
        'all_genes': collections.OrderedDict,
        'imgt_aligned_status': bool}
        attr_to_type.update({'cdr3_%s_aa_smat' % loci: parasail.bindings_v2.Matrix for loci in 'abgd'})
        attr_to_type.update({'cdr2_%s_aa_smat' % loci: parasail.bindings_v2.Matrix for loci in 'abgd'})
        attr_to_type.update({'cdr1_%s_aa_smat' % loci: parasail.bindings_v2.Matrix for loci in 'abgd'})
        attr_to_type.update({'pmhc_%s_aa_smat' % loci: parasail.bindings_v2.Matrix for loci in 'abgd'})
        attr_to_type.update({'cdr3_%s_aa_pw' % loci: np.ndarray for loci in 'abgd'})
        attr_to_type.update({'cdr2_%s_aa_pw' % loci: np.ndarray for loci in 'abgd'})
        attr_to_type.update({'cdr1_%s_aa_pw' % loci: np.ndarray for loci in 'abgd'})
        attr_to_type.update({'pmhc_%s_aa_pw' % loci: np.ndarray for loci in 'abgd'})
        attr_to_type.update({'dist_%s' % loci: pd.core.frame.DataFrame for loci in 'abgd'})


        with pd.HDFStore(hdf5_file) as store:

            # pull the metadata from cell_df
            metadata = store.get_storer('cell_df').attrs.metadata
            metadata = json.loads(metadata)

            for attr in metadata['attributes']:
                try:
                    if verbose:
                        print("SETTING {} AS {}".format(attr, attr_to_type[attr] ))
                except KeyError:
                    print("YOU TRIED TO SET {} BUT IT IS NOT A RECOGNIZED ATTRRIBUTE OF TCRrep".format(attr))
                    continue
                if not attr in store and not attr in metadata:
                    setattr(self, attr, None)
                    continue

                if attr_to_type[attr] is pd.core.frame.DataFrame:
                    setattr(self, attr, store[attr])
                elif attr_to_type[attr] is np.ndarray:
                    if 'pw' in attr:
                        tmp = store[attr].values
                        incoming = scipy.spatial.distance.squareform(np.squeeze(tmp),
                                                                     force='tomatrix',
                                                                     checks=False)
                        setattr(self, attr, incoming)
                    else:
                        setattr(self, attr, store[attr].values)
                elif attr_to_type[attr] in [int, str, list, bool, dict]:
                    setattr(self, attr, metadata[attr])
                

def load_hdf5(filename):
    """
    load a TCRrep saved to HDF5 using `save_as_hdf5`

    Parameters
    ----------
    filename : string
        Path to a HDF5 file containing a TCRrep

    Returns
    -------
    incoming : TCRrep
        a new TCRrep object"""

    """Instantiate an empty TCRrep and rebuild from file"""
    incoming = TCRrep(pd.DataFrame())
    incoming.rebuild_from_hdf5(filename)
    return incoming

def _map_gene_to_reference_seq(organism = "human",
                               gene= 'TRAV1-1*02',
                               cdr = 1,
                               ref = pb_cdrs):
    """
    get cdr amino acid seq from TCR V-gene name.


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
    """
    Use index_cols to group by and group identical entries. The input DataFrame
    must have a column 'count'.
    """
    clones = cell_df.groupby(index_cols)['count'].agg(np.sum).reset_index()
    return clones

def _compute_pairwise(sequences, metric='nw', processes=2, user_function=None, **kwargs):
    """Wrapper for pairwise.apply_pw_distance_metric_w_multiprocessing()

    Parameters
    ----------
    sequences : list
    metric : string
    processes : int
    user_function : function
    **kwargs : keyword arguments passed to metric function

    Returns
    -------
    pw_full_np : np.ndarray
        matrix of pairwise comparisons"""
    # processes = 1
    
    if metric == 'nw':
        metric_func = pwseqdist.metrics.nw_metric
    elif metric == 'hamming':
        metric_func = pwseqdist.metrics.nw_hamming_metric
    elif metric == 'tcrdist_cdr3':
        metric_func = pairwise.tcrdist_cdr3_metric
    elif metric in ['tcrdist_cdr1', 'tcrdist_cdr2', 'tcrdist_cdr2.5', 'tcrdist_pmhc']:
        metric_func = pairwise.tcrdist_cdr1_metric
    elif metric == 'hamming_aligned':
        metric_func = pwseqdist.metrics.hamming_distance
    elif not user_function is None:
        metric_func = user_function
    else:
        msg = 'repertoire._compute_pairwise: metric %s is not supported'
        raise ValueError(msg % metric)
    
    if metric in ['nw', 'hamming']:
        if not 'matrix' in kwargs:
            kwargs['matrix'] = 'blosum62'

    dvec = pwseqdist.apply_pairwise_sq(sequences.values, metric_func,
                                  ncpus=processes, **kwargs)
    """This may not be neccessary in all use cases of the distances.
    Consider returning the vector form."""
    pw_full_np = scipy.spatial.distance.squareform(dvec)
    """
    pw = pairwise.apply_pw_distance_metric_w_multiprocessing(
                    sequences = sequences, #! BUG FIX (sequences changed to unique_seqs)
                    metric = metric,
                    user_function = user_function,
                    processes= 1,
                    **kwargs)
    pw = pairwise._pack_matrix(pw)
    if not np.all(pw_full_np == pw):
        print('dvec', dvec)
        print(pw_full_np)
        print(pw)
        print(pw_full_np == pw)
        print('pwsd_outer', pw_full_np[1, 2])
        print('tcrdist_outer', pw[1, 2])
        print(sequences)
        print(sequences[1])
        print(sequences[2])
        print(kwargs['matrix'].name)
        print(metric)
    # pw_df = pd.DataFrame(pw, index=sequences, columns=sequences)

    #print(sequences.shape)
    #print(dvec.shape)
    #print(pw_full_np.shape)
    #print(dvec)
    
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
    pw_full_np = pw_full.values"""
    return pw_full_np
    
def _map_clone_df_to_TCRMotif_clone_df(df):
    """
    TODO: TEST REPLACEMENT WITH MUCH SIMPLER GENERIC
    mappers.generic_pandas_mapper(tr.clone_df, mappers.TCRrep_clone_df_to_TCRMotif_clone_df)

    THEN REMOVE AND REMOVE TESTS FROM test_repertoire_unit.py

    Converts clone_df DataFrame used in tcrdist2 to the input clones_df
    DataFrame required by TCRMotif().

    Parameters
    ----------
    df : DataFrame
        must contain columns ['subject','epitope', 'v_a_gene', 'j_a_gene',
                              'v_b_gene', 'j_b_gene', 'cdr3_a_aa','cdr3_b_aa']

    Returns
    -------
    df : DataFrame
        modified DataFrame with subset and renamed columns

    Example
    -------
    >>> df = pd.DataFrame([[1,2,3,4,5,6,7,8,9,10]],
                        columns = [ 'subject','epitope',
                                    'v_a_gene', 'j_a_gene',
                                    'v_b_gene', 'j_b_gene',
                                    'cdr3_a_aa','cdr3_b_aa',
                                    'cdr2_a_aa', 'cdr2_b_aa'])
    >>> print(_map_clone_df_to_TCRMotif_clone_df(df))
        subject  epitope  va_rep  ja_rep  vb_rep  jb_rep  cdr3a  cdr3b
    0        1        2       3       4       5       6      7      8
    """

    columns_conversion_dict =   {'subject'  : 'subject',
                                 'epitope'  : 'epitope',
                                 'v_a_gene' : 'va_rep',
                                 'j_a_gene' : 'ja_rep',
                                 'v_b_gene' : 'vb_rep',
                                 'j_b_gene' : 'jb_rep',
                                 'cdr3_a_aa': 'cdr3a',
                                 'cdr3_b_aa': 'cdr3b' }

    if not np.all([n in df.columns for n in columns_conversion_dict.keys()]):
        missing = [n for n in columns_conversion_dict.keys() if n not in df.columns ]
        raise KeyError("clone_df must have columns names: {}".format(" , ".join(missing)))

    df = df[columns_conversion_dict.keys()]\
        .rename(columns = columns_conversion_dict).copy()\

    return(df)


"""
Private Methods of TCRrep

Extended Summary
----------------
_validate_chains()
    raise ValueError is chains arg is mis-specified
_validate_chain()
    raise ValueError if chain arg is mis-specified
_validate_()
    raise TypeError if  is not pd.DataFrame
_initialize_chain_specific_attributes(self, chain)
    create chain specific attribute including setting default sub matrix
_get_smat()
    return smat given chain (e.g. "alpha") and index_col (e.g. "cdr2_a_aa")
_assign_pw_result()
    assign pw distance given chain (e.g. "alpha") and index_col (e.g. "cdr2_a_aa")


"""
