# base python imports (e.g., os, sys)
import os
import sys
# scientific python imports (e.g., numpy, pandas,)
import pandas as pd
# tcrdist specific imports
from .paths import path_to_current_db_files
from .all_genes import all_genes
from . import util
from . import find_cdr3_motifs_in_tcrdist2
from . import setup_db

class TCRMotif():
    """
    A class for motif finding, given a TCR clones file.

    Methods
    -------
    find_cdr3_motifs : function
        wrapper for tcrdist.find_cdr3_motifs_in_tcrdist2.find_cdr3_motif().
        It passes self.all_tcrs and self.ng_tcrs dictionaries along with
        user-supplied parameters. Results are assigned to the
        DataFrame self.motifs_df.
    generate_all_tcrs : function
        called on init, produces self.all_tcrs dict from self.clone_df
    generate_ng_tcrs : function
        called on init, produces self.ng_tcrs dict from the next generation sequence files.
        (these files must be previosly installed. (see tcrdist.setup_db)

    Attributes
    ----------

        User-Supplied Attributes
        ========================

    clones_df : DataFrame
        must contain the following minimum set of headers
        ['epitope','va_rep','ja_rep','vb_rep','jb_rep','cdr3a','cdr3b']
    organism  : string
        specifies the reference organism. e.g., "mouse", "human"
    chains : list
        list of strings specifying the TCR chains. e.g., ["A","B"]
    epitopes
        list of strings specifying the epitope(s) to consider e.g., ["PA"]
    motifs_df : DataFrame
        discovered motifs after successfully running self.find_cdr3_motifs()


        Default Attributes (Modification intended only for advanced users)
        ==================================================================

    min_count : int = 10
        see extend_motif()
    max_ng_lines : int
        num lines in the next-gen reference file to consider.
    max_motif_len : int = 100
        see extend_motif()
    nsamples :int = 25
        see extend_motif()
    min_expected : float= .25
        see extend_motif()
    max_overlap : float = 0.8
        see extend_motif()
    max_useful_expected : float= 1.0
        see extend_motif()
    chi_squared_threshold : float
        see extend_motif()
    chi_squared_threshold_for_seed : float
        see extend_motif()
    min_chi_squared_increase : float
        see extend_motif()
    min_extended_count_ratio : float
        see extend_motif()
    use_fake_seqs : bool
        see extend_motif()
    big : bool
        tcrdist1 uses True
    constant_seed : bool
        tcrdist1 uses True
    force_random_len  : bool
        tcrdist1 uses True
    verbose : bool
        outputs intermediate step to std.out
    nofilter  : bool
        tcrdist1 ?
    very_verbose  : bool
        outputs intermediate print statements
    test_random : bool
          tcrdist1 ?
    hacking : bool
          tcrdist1 ?
    unmask  : bool
          tcrdist1 ?
    use_groups  : bool
          tcrdist1 uses True
    begin : string
          designates beginning of string pattern, probably used in regex. e.g.,'^'
    end :
          designates beginning of end pattern, probably used in regex. e.g.,'$'
    X : string
          designates any capital letter, , probably used in regex. e.g.,'[A-Z]'
    dot : string
          designates a dot pattern, probably used in regex. e.g., '.'
    pseudocount : float
          e.g., 0.0


    Notes
    -----
    Initialization sets default parameter values, which can be manually
    overridden if desired.

    The default values, stored as TCRMotif attributes are converted to a
    dictionary and stored as self.params.

    This dictionary is important because it is passed
    as the **kwargs (keyword arguments) to TCRMotif.find_cdr3_motif().

    These settings control the behavior of the function
        trcdist.find_cdr3_motifs_in_tcrdist2.find_cdr3_motif()

    For instance, if the user wanted to change the chi_squared_threshold
    argument they would do so as follows:

    >>> motif = TCRMotif(clones_df = clones_df)
    >>> motif.chi_squared_threshold
    50.0
    >>> motif.chi_squared_threshold = 75.0
    >>> motif.chi_squared_threshold
    75.0
    """
    def __init__(self, clones_df, organism, chains, epitopes):
        """
        Parameters
        ----------
        clones_df : DataFrame
            must contain the following minimum set of headers
            ['epitope','va_rep','ja_rep','vb_rep','jb_rep','cdr3a','cdr3b']
        organism  : string
            specifies the reference organism. e.g., "mouse", "human"
        chains : list
            list of strings specifying the TCR chains. e.g., ["A","B"]
        epitopes
            list of strings specifying the epitope(s) to consider e.g., ["PA"]
        motifs_df : DataFrame
            discovered motifs after successfully running self.find_cdr3_motifs()
        """

        self.clones_df                = clones_df # provided by the dear-user
        #self.chain                    = #"A"
        self.chains                   = chains    #['A',"B"]
        self.organism                 = organism # "mouse"
        self.epitopes                 = epitopes  #['PA']

        self.all_tcrs                 = None  # set by self.generate_all_tcrs()
        self.ng_tcrs                  = None  # set by self.generate_ng_tcrs()
        self.motifs_df                = None  # generated by find_cdr3_motif

        # Default Control Attributes
        self.min_count                = 10
        self.max_ng_lines             = 10000
        self.max_motif_len            = 100
        self.nsamples                 = 25
        self.min_expected             = .25
        self.max_overlap              = 0.8
        self.max_useful_expected      = 1.0
        self.chi_squared_threshold    = 50.0
        self.chi_squared_threshold_for_seeds = 20.0
        self.min_chi_squared_increase = 25.0
        self.min_extended_count_ratio = 0.6
        self.use_fake_seqs            = False
        self.big                      = True
        self.constant_seed            = False
        self.force_random_len         = False
        self.verbose                  = True
        self.nofilter                 = False
        self.very_verbose             = False
        self.test_random              = False
        self.hacking                  = False
        self.unmask                   = True
        self.use_groups               = True
        self.begin                    = '^'
        self.end                      = '$'
        self.X                        = '[A-Z]'
        self.dot                      = '.'
        self.pseudocount              = 0.0

        # set params to param dict which will be used as **kwargs
        self.params= {}
        self._set_params()
        self.generate_all_tcrs() # This is the entire clones file
        self.generate_ng_tcrs()  # This is based on next gen sequence files

    def find_cdr3_motifs(self,
                         epitopes = None,
                         organism = None,
                         chains   = None,
                         all_tcrs = None,         # set by self.generate_all_tcrs()
                         ng_tcrs = None):

        if epitopes is None:
            epitopes = self.epitopes
        if organism is None:
            organism = self.organism
        if chains is None:
            chains = self.chains
        if all_tcrs is None:
            all_tcrs = self.all_tcrs
        if ng_tcrs is None:
            ng_tcrs  = self.ng_tcrs

        self._set_params() # makes sure that the latest params are used

        motifs_df = find_cdr3_motifs_in_tcrdist2.find_cdr3_motif(epitopes = epitopes,
                                                                organism = organism,
                                                                chains   = chains,
                                                                all_tcrs = all_tcrs ,         # set by self.generate_all_tcrs()
                                                                ng_tcrs  = ng_tcrs,
                                                                **self.params)
        if self.motifs_df is None:
            # update object attribute
            self.motifs_df = motifs_df
        else:
            # if motifs_df attribute already exists add to it, don't overwrite.
            self.motifs_df = self.motifs_df.append(motifs_df)
        return(motifs_df)

    def _set_params(self):
        """
        Sets self.params to a dictionary of these select parameters
        that are used throughout the class methods.

        If defaults are changed, TCRMotif._set_params() must be called again.

        It is called in find_cdr3_motif to ensure the most up to date settings
        are passed to

        """
        param_keys = ["min_count", "max_ng_lines", "max_motif_len",
                        "nsamples", "min_expected", "max_overlap",
                        "max_useful_expected", "chi_squared_threshold",
                        "chi_squared_threshold_for_seeds",
                        "min_chi_squared_increase",
                        "min_extended_count_ratio", "use_fake_seqs", "big",
                        "constant_seed", "force_random_len", "verbose",
                        "nofilter", "very_verbose", "test_random",
                        "hacking", "unmask", "use_groups", "begin",
                        "end", "X", "dot", "pseudocount"]
        full_params = self.__dict__

        self.params =  {k:full_params[k] for k in param_keys}
        return {k:full_params[k] for k in param_keys}

    def generate_all_tcrs(self, clones_df = None, organism = None):
        """
        Parameters
        ----------
        clones_df : DataFrame
        organism : string

        Returns:
        all_tcrs : dict
        """
        if clones_df is None:
            clones_df = self.clones_df
        if organism is None:
            organism = self.organism

        df = self._populate_all_tcrs_df_from_clones_df(clones_df = clones_df, organism = organism)
        all_tcrs = self._convert_all_tcrs_df_to_all_tcrs_dict(df)
        self.all_tcrs = all_tcrs
        return(all_tcrs)


    def generate_ng_tcrs(self,
                         chains = None,
                         organism = None,
                         max_ng_lines = None):
        """
        populate next gen tcrs dict, using next generation sequence data

        Parameters
        ----------
        chains: list
            list of strings default is ["A","B"]
        organism : string
            either mouse of human
        max_ng_lines : int
            number of lines in next-gen file to consider

        Returns
        -------
        ng_tcrs : dict
            dict (chain) of dict(v_gene) of dicts (j_gene) pointing to list of tuples

        Raises
        ------
        OSError if required next gen reference files are not in tcrdist/db/alphabeta_db/

        Example
        -------
        ng_tcrs = populate_nextgen_tcrs_dict()
        ng_tcrs['A']['TRAV14-2*01']['TRAJ12*01']
        [('CAASVGGYKVVF', 'tgtgcagcaagtgttggaggctataaagtggtcttt'),
         ('CAASPGTGGYKVVF', 'tgtgcagcaagtcccgggactggaggctataaagtggtcttt'),
         ....
         ]
        """
        if chains is None:
            chains = self.chains
        if organism is None:
            organism = self.organism
        if max_ng_lines is None:
            max_ng_lines = self.max_ng_lines

        ng_tcrs = {k:{} for k in chains}

        for ab in chains:
            ng_logfile = os.path.join(path_to_current_db_files(), 'new_nextgen_chains_{}_{}.tsv'.format(organism,ab)) # - INPUT FROM df files

            if not os.path.isfile(ng_logfile):
                raise OSError('find_cdr3_motifs.py: missing next-gen chains file {}'.format(ng_logfile))

            counter=0
            num_chains=0
            ab_chains = {}

            for line in open(ng_logfile,'r'):
                counter+=1
                l = line[:-1].split('\t')
                if counter==1:
                    assert l==['v_reps','j_reps','cdr3','cdr3_nucseq']
                    continue
                #if not counter%1000000:#Log(`counter`+' '+`num_chains`+' '+ng_logfile)
                if max_ng_lines and counter > max_ng_lines:
                    break
                v_reps = set( ( util.get_mm1_rep(x,organism) for x in l[0].split(',') ) )
                j_reps = l[1].split(',')
                cdr3, cdr3_nucseq = l[2:4]

                ## now add to the different places
                for v_rep in v_reps:
                    for j_rep in j_reps:
                        if v_rep not in ab_chains: ab_chains[v_rep] = {}
                        if j_rep not in ab_chains[v_rep]: ab_chains[v_rep][j_rep] = []
                        ab_chains[v_rep][j_rep].append( (cdr3, cdr3_nucseq ))


                num_chains += 1

                ng_tcrs[ab] = ab_chains
        self.ng_tcrs = ng_tcrs
        return ng_tcrs


    def _populate_all_tcrs_df_from_clones_df(self, clones_df, organism):
        """
        Parameters
        ----------
        clones_df : DataFrame
        organism  : string

        Returns
        -------
        all_tcr : DataFrame

        """
        all_tcrs_df = clones_df[['epitope','va_rep','ja_rep','vb_rep','jb_rep','cdr3a','cdr3b']]
        all_tcrs_df = all_tcrs_df.rename(columns = {'va_rep':'va','ja_rep':'ja','vb_rep':'vb','jb_rep':'jb'})
        all_tcrs_df['va_rep'] = all_tcrs_df['va'].apply(lambda k: all_genes[organism][k].mm1_rep).copy()
        all_tcrs_df['ja_rep'] = all_tcrs_df['ja'].apply(lambda k: all_genes[organism][k].mm1_rep).copy()
        all_tcrs_df['vb_rep'] = all_tcrs_df['vb'].apply(lambda k: all_genes[organism][k].mm1_rep).copy()
        all_tcrs_df['jb_rep'] = all_tcrs_df['jb'].apply(lambda k: all_genes[organism][k].mm1_rep).copy()
        return(all_tcrs_df)

    def _convert_all_tcrs_df_to_all_tcrs_dict(self, df, tup_order = ["va", "ja",
                                                           "vb", "jb",
                                                           "cdr3a", "cdr3b",
                                                           "cdr3a", "cdr3b",
                                                           "va_rep", "ja_rep",
                                                           "vb_rep", "jb_rep"]):
        """
        Converts a panda DataFrame to dictionary of tuples keyed on epitope.
        Tuple order is specified in tup_order

        Parameters
        ----------
        df : DataFrame
        tup_order : list


        Returns
        -------
        all_tcrs : dict

        Notes
        -----
        This directly mimics the pre-existing function that returns a tuple
        ( va, ja, vb, jb, cdr3a, cdr3b, cdr3a, cdr3b, va_rep, ja_rep, vb_rep, jb_rep )

        """
        all_tcrs = {}
        for i, r, in df.iterrows():
            ep = r['epitope']
            row_values = r[tup_order].to_list()
            row_tupple = tuple(row_values)
            if ep not in all_tcrs.keys():
                all_tcrs[ep] = []
            all_tcrs[ep].append(row_tupple )
        return(all_tcrs)

    def validate_ng_files_installed(self, organism = None, chains = None):
        """
        WHEN IMPLEMENTED THIS WILL HAVE TO BE INTEGRATED INTO THE INIT METHOD.
        TODO: MAKE THE DOWNLOAD SCRIPT RECURSIVE AND MAKE IT WORK FOR WINDOWS,
        OR OFFER DOWNLOAD INSTRUCTIONS FOR WINDOWS USERS.
        MAKE SURE THIS WOULD WORK FOR A USER WHO PIP INSTALLED THIS INSTEAD OF
        CLONED IT. THIS IS A REALLY GOOD EXAMPLE OF WHERE DOCKER WILL BE THE
        WAY FORWARD.
        """
        if chains is None:
            chains = self.chains
        if organism is None:
            organism = self.organism

        for ab in chains:
            ng_logfile = os.path.join(path_to_current_db_files(), 'new_nextgen_chains_{}_{}.tsv'.format(organism,ab)) # - INPUT FROM df files
            if not os.path.isfile(ng_logfile):
                #raise OSError('WARNING:: find_cdr3_motifs.py: missing next-gen chains file {}'.format(ng_logfile))
                print("Do you want to install new_nextgen_chains_{}_{}.tsv?[Y/N]".format(organism,ab))
                user_response = input()
                if user_response.upper() not in ["Y","N", "YES","NO"]:
                    print("You must supply 'Y' or 'N'")
                    user_response = input()
                elif user_response.upper() in ["Y","YES"]:
                    setup_db.install_nextgen_data_to_db(download_file = 'new_nextgen_chains_{}_{}.tsv'.format(organism,ab))
                else:
                    raise OSError('find_cdr3_motifs.py: missing next-gen chains file {}'.format(ng_logfile))
