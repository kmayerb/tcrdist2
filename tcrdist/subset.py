# base python imports
import re
import random

# python env imports
import numpy as np
import pandas as pd

# tcrdist2 imports
from .all_genes import all_genes
from . import amino_acids
from . import basic
from . import logo_tools
from . import mappers
from . import paths
from . import random
from . import rmf
from . import svg_basic
from . import tcr_sampler
from . import util
from collections import namedtuple
from .storage import StoreIO
from .storage import StoreIOMotif
from .storage import StoreIOEntropy
from .cdr3_motif import TCRMotif

#from tcrdist import parse_tsv

class TCRsubset():
    """
    A class dedictated to repertoire subset analysis,
    particularly relative entropy motifs

    """
    def __init__(self,
                 clone_df,
                 organism = None,
                 epitopes = None,
                 epitope  = None,
                 chains   = None,
                 nbr_dist = None,
                 dist_a   = None,
                 dist_b   = None,
                 dist_g   = None,
                 dist_d   = None):

        # init param attributes
        self.epitopes = epitopes
        self.epitope  = epitope
        self.clone_df = clone_df
        self.organism = organism
        self.chains   = chains
        self.dist_a   = dist_a
        self.dist_b   = dist_b
        self.dist_g   = dist_g
        self.dist_d   = dist_d
        self.nbr_dist = nbr_dist # !!!!!

        # apply defaults, when None are supplied at initiation
        if self.nbr_dist is None:
            self.nbr_dist = 100.0

        # placeholders for internally generated attributes
        self.all_tcrs      = None
        self.ng_tcrs       = None
        self.all_rep2label_rep = None
        self.all_rep2label_rep_color = None
        self.motif_df = None

        # Epitope Specific
        self.tcrs          = None
        self.rep2label_rep = None
        self.rep2label_rep_color = None
        self.motifs        = None

        # Motif Specific
        self.all_neighbors = None
        self.vl     = None
        self.jl     = None
        self.vl_nbr = None
        self.jl_nbr = None


        # shared class resources for validation
        self.chain_to_dist = {"A": "dist_a",
                              "B": "dist_b",
                              "G": "dist_g",
                              "D" :"dist_d"}

        # Validation
        self._validate_chains()
        self._validate_organism()

        # Generation
        self._generate_tcrs(self.epitope)
        self._generate_all_neighbors()
        self._generate_ng_tcrs()

        # Validation
        self._validate_all_dists_and_tcrs_match()

    def tcr_motif_clones_df(self):
        """
        Use this function to create a clones_df input appropriate to TCRMotif.

        It make use of a mapper to ensure proper columns and column names

        Example
        -------
        TCRMotif(clones_df = TCRSubset.tcr_motif_clones_df())
        """
        return mappers.generic_pandas_mapper(self.clone_df,
                                             mappers.TCRsubset_clone_df_to_TCRMotif_clone_df)

    def find_motif(self):
        """
        Create a TCRMotif_instance using subset organism, chains, and epitopes.
        runs TCRMotif.find_motif. Warning this can take 5-10 minutes per chain.

        Returns
        -------
        motif_df : DataFrame
        """
        TCRMotif_instance = TCRMotif( self.tcr_motif_clones_df(),
                                      organism = self.organism,
                                      chains = self.chains,
                                      epitopes = self.epitopes)
        print("SEARCHING FOR MOTIFS, THIS CAN TAKE 5-10 minutes")
        TCRMotif_instance.find_cdr3_motifs()
        motif_df = TCRMotif_instance.motif_df.copy()
        self.motif_df = motif_df
        return motif_df

    def loop_through_motifs(self):
        """
        NOT READY BUT WILL USE TO LOOP THROUGH
        """
        if motif_df is None:
            motif_df = self.motif_df
        if motif_df is None:
            raise ValueError("TCRsubest.motif_df DataFrame is empty. Load one or try TCRsubest.motif_df.find_motif()")
        pass
        svg_list = []
        for i,row in motif_df.iterrows():
            StoreIOMotif_instance = StoreIOMotif(**row)
            self.analyze_motif(s = StoreIOMotif_instance)
            self.analyze_matches(s = StoreIOMotif_instance)
            svg = plot_pwm(StoreIOMotif_instance, create_file = False, my_height = 200, my_width = 600)
            svg_list.append(svg)
        return svg_list

    def eval_motif(self, row):
        """
        eval motif wraps functions for evaluating a row of the motif_df DataFrame

        row : OrderedDict
            from a row of motif_df
            i = 0; row = tm.motif_df.iloc[i,:].to_dict()

        Returns
        -------
        StoreIOMotif_instance : StoreIOMotif

        Raises
        ------
        TypeError
            if row variables cannot be coerced to correct types
            (see: StoreIOMotif_instance._coerce_attrs() )

        Notes
        -----
        The steps shown above are:

        1. Initialize instance of the information carrier class StoreIOMotif
        2. Analyze_motif (ts.analyze_motif) to determine matches and neighbors
        of matches, with new attributes are appended to the StoreIOMotifinstance.
        3. Analyze matches (ts.analyze_matches) to identify the the relative
        entropy between position wise matrices
        """
        # 1
        StoreIOMotif_instance = StoreIOMotif(**row)
        StoreIOMotif_instance._coerce_attrs()
        assert StoreIOMotif_instance._validate_attrs()
        # 2
        self.analyze_motif(s = StoreIOMotif_instance)
        # 3
        self.analyze_matches(s = StoreIOMotif_instance)
        return StoreIOMotif_instance


    def _load_motifs_from_file(self):
        #motif_fn = 'mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones_cdr3_motifs_PA.log'
        #motif_fh =open(motif_fn, 'r')
        pass


    def _load_motifs_from_dataframe(self):
        pass


    def analyze_motif(self, s, tcrs = None, all_neighbors = None):
        """
        Parameters
        ----------
        s : tcrdist.storage.StoreIOMotif
            StoreIOMotif instance
        tcrs : dict
            if omitted, default is to self.tcrs
        all_neighbors : dict
            if omitted, default is to self.all_neighbors
        Returns
        -------
        s : tcrdist.storage.StoreIOMotif
            StoreIOMotif instance with new attributes added
        Notes
        -----
        The primary function of this script is to add the following attributes
        to a StoreIOMotif instance:
            s.showmotif
            s.vl_nbr
            s.jl_nbr
            s.vl
            s.jl
            s.matches - tcrs that perfectly match regex
            s.nbr_matches - tcr hat perfectly match regex + neigbors
                (tcrdist<self.nbr_dist)
            s.matched_tcrs_plus_nbrs - index positions
            s.matched_tcrs - index positions

        Developers Notes - Opportunity for Modularity

        A goal of refactoring is to avoid class bloat.

        * All self attributes are defined at the top of this function.
        * Nothing is assigned back to self.

        Therefore, it is possible to refactor, by removing
        the body of this function from the
        class and call it from a clean wrapper.

        !! self._get_counts_lists_from_tcr_indices is called and would also need
        to be separated from the class.
        """
        # Pull Reference objects. These are in variant of the motif considered
        if tcrs is None:
            tcrs = self.tcrs

        num_tcrs = len(tcrs)

        if all_neighbors is None:
            all_neighbors = self.all_neighbors

        # Unpack motif specific variables. These depend on the motif being
        # considered. String variables are unpacked from s -
        # the storage_io_motif_object, and the types are enforced.

        count          = int(s.count)
        expect_random  = float(s.expect_random)
        expect_nextgen = float(s.expect_nextgen)
        chi_squared    = float(s.chi_squared)
        nfixed         = int(s.nfixed)
        showmotif      = str(s.showmotif)
        num            = int(s.num)
        othernum       = int(s.othernum)
        overlap        = int(s.overlap)
        ep             = str(s.ep)
        ab             = str(s.ab)
        nseqs          = int(s.nseqs)
        v_rep_counts   = str(s.v_rep_counts)
        j_rep_counts   = str(s.j_rep_counts)

        showmotif = list(showmotif)

        expected_fraction = max(expect_random, expect_nextgen) / num_tcrs
        motif = [amino_acids.groups[x] for x in showmotif]

        # motif e.g.,^.ALGaGaN as re pattern e.g., ^[A-Z]ALG[AGSP]G[AGSP]N
        prog = re.compile(''.join(motif))

        total = 0
        matches = []                # Comprised of prefect matches to re.motif
        nbr_matches = []            # Comprised of perfect matches and neighbors of perfect matches
        matched_tcrs = []           # Comprised of index positions of perfect matched
        matched_tcrs_plus_nbrs = [] # Above plus, index positions of neighbors of perfect matches

        # < tcrs : list, tuple > are all epitope-specific tcrs, we search them all
        for ii, tcr in enumerate( tcrs ):
            # select the cdr3 (e.g., CAMRGNSGGSNYKLTF) and cdr3_nucseq_src ('V', 'V', ..., 'J, 'J')
            # TODO : GENERALIZE THIS SO THAT IT DOES NOT SEARCH FOR "A"
            if ab in ["A","G"]:
                cdr3,cdr3_nucseq_src = tcr[4], tcr[10]
            else:
                cdr3,cdr3_nucseq_src = tcr[5], tcr[11]

            # search for the motif re in each cdr3
            m = prog.search(cdr3)
            # if found
            if m:
            # < mseq : str > portion of the cdr3 amino acid that matches regex pattern
                mseq = cdr3[ m.start():m.end() ]
            # < nseq : str > source (V, N, or J ) of the cdr3 nucleotide that matches the regex pattern
                nseq_src = cdr3_nucseq_src[ 3*m.start():3*m.end() ]
            # < positions : range >
                positions = range(m.start(),m.end())
            # < rpositions > reverse position relative to the cdr3 end
                rpositions = [len(cdr3)-1-x for x in positions]
            # < matches : list > outside the loop contains mseq,nseq,positions,rpositions
                matches.append( (mseq,nseq_src,positions,rpositions) )
            # < matched_tcrs: list > outside the loop contains the index number of the matching tcr
                matched_tcrs.append( ii )

            # < all_neighbors > in same order as tcrs
                for nbr in all_neighbors[ab][ii]:
                    if nbr not in matched_tcrs_plus_nbrs:
                    # < nbr_tcr > neighbor to a perfect match pulled from < tcr >
                        nbr_tcr = tcrs[nbr]
                        if ab in ['A',"G"]:
                            nbr_cdr3, nbr_cdr3_nucseq_src = nbr_tcr[4], nbr_tcr[10]
                        else:
                            nbr_cdr3, nbr_cdr3_nucseq_src = nbr_tcr[5], nbr_tcr[11]
                    # !! only if nbr and principal cdr3 are same length will they be included
                        if len(nbr_cdr3) == len(cdr3):
                            matched_tcrs_plus_nbrs.append( nbr )
                            nbr_mseq = nbr_cdr3[ m.start():m.end() ]
                            nbr_nseq = nbr_cdr3_nucseq_src[ 3*m.start():3*m.end() ]
                            nbr_matches.append( (nbr_mseq,nbr_nseq,positions,rpositions) )


                total += 1

        vl_nbr, jl_nbr = self._get_counts_lists_from_tcr_indices(matched_tcrs_plus_nbrs, ab)
        vl, jl = self._get_counts_lists_from_tcr_indices(matched_tcrs, ab)

        # send outputs back to input object IO
        s.showmotif              = showmotif # now returned as a list
        s.vl_nbr                 = vl_nbr
        s.jl_nbr                 = jl_nbr
        s.vl                     = vl
        s.jl                     = jl
        s.matches                = matches
        s.nbr_matches            = nbr_matches
        s.matched_tcrs_plus_nbrs = matched_tcrs_plus_nbrs
        s.matched_tcrs           = matched_tcrs
        return(s)


    def analyze_matches(self, s):
        """
        Parameters
        ----------
        s : StorageIOMotif
            a StorageIOMotif without entropy information

        Returns
        -------
        s : StorageIOMotif
            Updated StorageIOMotif instance w/ entropy attribute
            which points to a StorageIOEntropy instance
        """
        StorageIOEntropy_instance = self._analyze_matches_using_ngseqs(
                                        matches = s.nbr_matches ,
                                        matched_tcrs = s.matched_tcrs_plus_nbrs,
                                        ab = s.ab,
                                        epitope = s.ep,
                                        showmotif = s.showmotif,
                                        tcrs = self.tcrs,
                                        ng_tcrs = self.ng_tcrs,
                                        num_nextgen_samples = 100,
                                        junction_bars = True,
                                        junction_bars_order = {'B': ['V','N1','D',
                                                                     'N2','J'],
                                                             'A': ['V','N','J'] },
                                        min_prob_for_relent_for_scaling = 1e-3,
                                        max_column_relent_for_scaling = 3.0)
        # Add storage_io_entropy_instance to storage_io_motif_instance
        s.entropy = StorageIOEntropy_instance
        return(s)


    def _generate_all_tcrs(self):
        """
        generates self.all_tcrs attribute

        Parameters
        ----------
        epitope : str

        See extensive documentation associated with:
        rmf._generate_tcrs_dict_from_clones_dataframe()

        """
        # implement backwards mapping to names recognized by tcrdist motif routine
        clone_df = mappers.generic_pandas_mapper(self.clone_df, mappers.tcrdist2_to_tcrdist_clone_df_mapping)

        # *_ dumps the 2nd and 3rd parts of the tuple,
        all_tcrs, all_rep2label_rep, all_rep2label_rep_color =\
            rmf._generate_tcrs_dict_from_clones_dataframe(clone_df,
                                                     epitopes = self.epitopes,
                                                     organism = "mouse",
                                                     return_as_tuple= True)
        self.all_tcrs                = all_tcrs
        self.all_rep2label_rep       = all_rep2label_rep
        self.all_rep2label_rep_color = all_rep2label_rep_color

        return all_tcrs

    def _generate_tcrs(self, epitope):
        """
        generates self.tcrs attribute

        Parameters
        ----------
        epitope : str

        See extensive documentation associated with:
        rmf._generate_tcrs_dict_from_clones_dataframe()

        """
        tcrs = self._generate_all_tcrs()[epitope]

        self.tcrs = tcrs
        # these were stored_when running _generate_all_tcrs
        self.rep2label_rep       = self.all_rep2label_rep[epitope]
        self.rep2label_rep_color = self.all_rep2label_rep_color[epitope]

        return tcrs

    def _generate_all_neighbors(self, nbr_dist = None):
        """
        generates self.all_neighbors attribute, using rmf.generate_all_nbr_from_dataframe

        Returns
        -------
        all_neighbors : dict
            dictionary keyed on chain (e.g., 'A', 'B'),
            containing list of lists for each position i in the list,
            the ith list contains the index position of the tcrs
            within (nbr_dist) distance from the ith tcr

        Assigns
        -------
        self.all_neighbors
            dictionary keyed on chain (e.g., 'A', 'B'),
            containing list of lists for each position i in the list,
            the ith list contains the index position of the tcrs
            within (nbr_dist) distance from the ith tcr

        Raises
        ------
        OSError if a required dist_x is not loaded in self

        AssertionError if length of all_nbr does not equal the
            row or column dimension of distance matrix.

        Notes
        -----
        rmf module refers to read_motif functions, re-factored from tcrdist1

        """
        if nbr_dist is None:
            nbr_dist = self.nbr_dist

        all_neighbors = {chain: None for chain in self.chains}

        for chain in self.chains:

            # converts chain letter 'X' to 'dist_x'
            dist_name = self.chain_to_dist[chain]

            # gets distance matrix matching chain
            if getattr(self, dist_name) is not None:
                dist_df = getattr(self, dist_name)
            else:
                raise OSError("{} not loaded".format(dist_name))

            # produces list of lists [[],[],...], length must equal dimensions of dist_df
            all_nbr = rmf.generate_all_nbr_from_dataframe(dist_df = dist_df, nbr_distance = nbr_dist)

            assert len(all_nbr) == dist_df.shape[0]
            assert len(all_nbr) == dist_df.shape[1]

            all_neighbors[chain] = all_nbr

        self.all_neighbors = all_neighbors
        return all_neighbors

    def _generate_ng_tcrs(self):
        self.ng_tcrs = rmf._generate_read_motif_ng_tcrs_dict(chains = self.chains)

    def _validate_chains(self):
        """
        Check that chains are a valid selection
        """
        valid_chains = ["A","B","G","D"]
        for chain in self.chains:
            if chain not in valid_chains:
                raise ValueError('chains must be one of ["A","B","G","D"]')

    def _validate_organism(self):
        """
        check that organism is valid string
        """
        valid_organisms = ["mouse", "human"]
        if self.organism not in valid_organisms:
            raise ValueError('organism must be one of ["mouse", "human"]')

    def _validate_all_dists_and_tcrs_match(self):
        """
        Checks that all chains in self.chains
        that the order of the distance matrices and
        tcrs list attribute match perfectly
        """

        for chain in self.chains:
            d = self.chain_to_dist[chain]
            self._validate_dist_and_tcrs_match(d)

    def _validate_dist_and_tcrs_match(self, dist):
        """
        Checks for a given chain distance DataFrame
        that the order of the distance matrix and
        tcrs list attribute match perfectly
        """
        if getattr(self, dist) is not None:
            x = getattr(self, dist)
        else:
            raise OSError("{} not loaded".format(dist))
        if not isinstance(x, pd.DataFrame):
            raise TypeError("distances must be DataFrames in TCRsubset")

        dist_order_row = x.index
        dist_order_col = x.columns

        tcrs_order = [x[-1]['clone_id'] for x in self.tcrs]

        if not np.all(dist_order_row == dist_order_col ):
            raise ValueError("index and columns must match for {}".format(dist))
        if not np.all(dist_order_row == tcrs_order ):
            raise ValueError("dist order and tcrs order must match")

    def _analyze_matches_using_ngseqs(self,
                                     matches,
                                     matched_tcrs,
                                     ab,
                                     epitope,
                                     showmotif,
                                     tcrs,
                                     ng_tcrs,
                                     num_nextgen_samples = 100,
                                     junction_bars = True,
                                     junction_bars_order = { 'B': ['V','N1','D','N2','J'], 'A': ['V','N','J'] },
                                     min_prob_for_relent_for_scaling = 1e-3,
                                     max_column_relent_for_scaling = 3.0):
        """
        Analyzes the motif matching sequences and compares them to reference nextgen
        sequences to discover (relative entropy) how the positive wise frequency
        distribtion of the motif is different from a reference probability
        distribution comprised of CDR3s from the same VJ-gene usage.

        Refer Questions about the Algorithm to Phil Bradley.

        Parameters
        ----------
        matches : list
            list containing motif matches [((mseq,nseq_src,positions,rpositions)), ]
                mseq : e.g., CALGGGSN
                nseq_src : e.g.,['V', 'V', 'V', 'V', ..., 'N', 'N', 'N', 'N', 'N', 'J', 'J', 'J', 'J', 'J', 'J']
                positions : range(0, 8),
                rpositions : [12, 11, 10, 9, 8, 7, 6, 5]
        matched_tcrs : list
            list [16, 49, 92, 112,...] index of those tcrs
        ab : str
            indicates chain
        epitope : str

        showmotif : list

        tcrs : list
            list of tcr information
        ng_tcrs : dict
            dict of form {chain:{v:{j:[]}}}
        num_nextgen_samples

        junction_bars : bool

        min_prob_for_relent_for_scaling : float

        max_column_relent_for_scaling : float

        Returns
        -------
        [dict, dict, dict, dict, dict, dict, dict, dict, list, list, int, int]
        pwm : dict
        npwm : dict
        ng_lenpwm : dict
        ng_fwdpwm : dict,
        ng_revpwm : dict
        fwdpwm : dict
        revpwm : dict
        scale_by_relent   : dict
        ng_fwdseq_reps,   : list
        ng_lenseq_reps    : list
        len( ng_lenseqs ) : int
        len( ng_fwdseqs ) : int

        Notes
        -----

        TODO: Finish a complete explanation of what the heck is going on here

        1. given a motif - ^.ALGaGaN (read in from from the motif_file)

        2. a regex motif - ^[A-Z]ALG[AGSP]G[AGSP]N is generated to be permissive in lowercase position)

        3. < matches >  contains information on those tcrs that were recognized by the regex motif

        4. < matched tcrs > contains index position of the those matches
        """
        ng_lenseqs = []
        ng_fwdseqs = []
        ng_revseqs = []

        ng_fwdseq_reps = []
        ng_lenseq_reps = []
        matched_reps = []

        seen = set() ## no repeats of ngseqs
        seen_samelen = set() ## no repeats of ngseqs

        for (mseq,nseq,positions,rpositions),ii in zip( matches, matched_tcrs ):
            tcr = tcrs[ii]
            if ab == 'A':
                my_cdr3,vrep,jrep = tcr[4:5]+tcr[6: 8]
            else:
                my_cdr3,vrep,jrep = tcr[5:6]+tcr[8:10]
            matched_reps.append( ( vrep, jrep ) )

            mylen = len(my_cdr3)
            if vrep in ng_tcrs[ab] and jrep in ng_tcrs[ab][vrep]:
                ngl = [ x for x in ng_tcrs[ab][vrep][jrep] if x not in seen ]
                if not ngl:
                    pass
                    #print('empty ngl!')

                for ngseq in random.sample( ngl, min(num_nextgen_samples,len(ngl)) ):
                    seen.add(ngseq)
                    (cdr3,cdr3_nucseq) = ngseq
                    L = len(cdr3)
                    fseq = ''
                    rseq = ''
                    for pos in positions:
                        if pos>=L:
                            fseq += '-'
                        else:
                            fseq += cdr3[pos]
                    for pos in rpositions:
                        if pos>=L:
                            rseq += '-'
                        else:
                            rseq += cdr3[L-1-pos]
                    ng_fwdseqs.append(fseq)
                    ng_revseqs.append(rseq)
                    ng_fwdseq_reps.append( ( vrep, jrep ) )

                ## cdr3s with the same length
                ngl_samelen = [ x for x in ng_tcrs[ab][vrep][jrep] if len(x[0]) == mylen and x not in seen_samelen ]
                if not ngl_samelen:
                    pass
                    #print('empty ngl_samelen!')
                for ngseq in random.sample( ngl_samelen, min(num_nextgen_samples,len(ngl_samelen))):
                    seen_samelen.add( ngseq )
                    cdr3 = ngseq[0]
                    ng_lenseqs.append( ''.join( [ cdr3[x] for x in positions ] ) )
                    ng_lenseq_reps.append( ( vrep, jrep ) )

        pwm = logo_tools.create_protein_pwm_from_sequences( [x[0] for x in matches ])

        npwm_alphabet = junction_bars_order[ab] if junction_bars else ['V','N','D','J']
        npwm = logo_tools.create_pwm_from_sequences( [x[1] for x in matches ], npwm_alphabet )

        #nbr_pwm = logo_tools.create_protein_pwm_from_sequences( [x[0] for x in nbr_matches ])
        #nbr_npwm = logo_tools.create_pwm_from_sequences( [x[1] for x in nbr_matches ], ['V','D','J','N'] )

        if ng_lenseqs:
            ng_lenpwm = self.create_wtd_pwm_from_sequences( ng_lenseqs, amino_acids.amino_acids+['-'], matched_reps, ng_lenseq_reps )
        else:
            ng_lenpwm = 0

        ng_fwdpwm = self.create_wtd_pwm_from_sequences( ng_fwdseqs, amino_acids.amino_acids+['-'], matched_reps, ng_fwdseq_reps )
        ng_revpwm = self.create_wtd_pwm_from_sequences( ng_revseqs, amino_acids.amino_acids+['-'], matched_reps, ng_fwdseq_reps )

        N = len(pwm)
        fwdpwm = {}
        revpwm = {}
        for i in range(N):
            fwdpwm[i] = {}
            revpwm[i] = {}
            incrememnt = 1.0/len(matches)
            for pos in [x[2][i] for x in matches]:
                fwdpwm[i]['pos'] = fwdpwm[i].get('pos',0)+incrememnt
            for pos in [x[3][i] for x in matches]:
                revpwm[i]['pos'] = revpwm[i].get('pos',0)+incrememnt

        ## look at relative entropies between nbrpwm and the fwd and rev pwms
        ## not nbr anymore since this is a subroutine
        ##
        scale_by_relent = {}
        for i in range(N):
            relents=[]
            for control_pwm in [ ng_fwdpwm[i], ng_revpwm[i] ]:
                relent = 0.0
                for a,pa in pwm[i].items():
                    if pa>= min_prob_for_relent_for_scaling:
                        qa = max(min_prob_for_relent_for_scaling, control_pwm.get(a,min_prob_for_relent_for_scaling))
                        paqa = np.log2(pa/qa)
                        relent += pa * paqa
                relents.append( relent )
            scale_by_relent[i] = max(0.,min(1., min(relents)/max_column_relent_for_scaling) )
            print('RE {:2d} {:5.2f} {:5.2f} {:5.2f} {} {} {}'.format( i, min(relents), relents[0], relents[1], ab, epitope, ''.join(showmotif) ))

        result_analyze_matches = (pwm, npwm, ng_lenpwm, ng_fwdpwm, ng_revpwm,\
                                  fwdpwm, revpwm, scale_by_relent, ng_fwdseq_reps,\
                                  ng_lenseq_reps, len( ng_lenseqs ), len( ng_fwdseqs))

        # result_analyze_matches is messy tuple containing
        # ([dict, dict, dict, dict, dict, dict, dict, dict, list, list, int, int])
        # So First, check that types are correct
        for i,t in enumerate([dict, dict, dict, dict, dict, dict, dict, dict, list, list, int, int]):
            assert isinstance(result_analyze_matches[i], t)
        field_names = ["pwm", "npwm", "ng_lenpwm", "ng_fwdpwm", "ng_revpwm", "fwdpwm",
              "revpwm", "scale_by_relent", "ng_fwdseq_reps", "ng_lenseq_reps",
              "num_ng_lenseqs", "num_ng_fwdseqs"]
        ES = namedtuple('ES', field_names)
        # Use namedtuple as an efficient way to drop tuple outputs
        # into a ordered dictionary
        result_analyze_matches_dict = ES(*result_analyze_matches)._asdict()
        # Now We Return Result in a Tidy Object: a StoreIOEntropy instance
        store_io_entropy = StoreIOEntropy(**result_analyze_matches_dict)
        return(store_io_entropy)


    def create_wtd_pwm_from_sequences(self, seqs, alphabet, target_reps, reps ):
                assert len(seqs) == len(reps)
                num_target_reps = len(target_reps)

                for bigrepeat in range(10):
                    reppair_wts = {}
                    if bigrepeat==0:
                        for rp in reps:
                            reppair_wts[rp] = 1.0 ## starting guess
                    else:
                        for rp in reps:
                            reppair_wts[rp] = 0.75 + 0.5 * random.random() ## starting guess

                    prev_dev = 1e6
                    for repeat in range(100):

                        ## what's the deviation
                        dev = 0.0
                        for ii in range(2):
                            ii_target_reps = [x[ii] for x in target_reps]
                            ii_reps = [x[ii] for x in reps]

                            scale_factor = float( len(reps ) )/ len(target_reps)

                            counts = {}
                            for rp in reps:
                                counts[rp[ii]] = counts.get(rp[ii],0) + reppair_wts[rp]

                            for rep,count in counts.items():
                                desired_count = scale_factor * ii_target_reps.count(rep)
                                dev += abs( desired_count - count )
                                fac = float(desired_count)/count
                                adjust = fac**(0.25)

                                #print 'desired_count:',desired_count,'count:',count,'fac:',fac,'adjust:',adjust,rep

                                for rp in reppair_wts:
                                    if rp[ii] == rep:
                                        reppair_wts[rp] *= adjust
                        #print 'repeat:',repeat,'dev:',dev
                        if abs(prev_dev-dev)<1e-3 and dev<1e-1:
                        #if abs(dev)<1e-1:
                            break
                        prev_dev = dev

                    #print 'final_dev:', bigrepeat,dev
                    if dev<1e-1:
                        break


                L = len(seqs[0])
                pwm = {}
                for i in range(L):
                    pwm[i] = dict(zip(alphabet,[0.0]*len(alphabet)))

                for seq,rp in zip( seqs, reps ):
                    assert len(seq) == L
                    seqwt = reppair_wts[ rp ]
                    #print seq, rp, seqwt
                    for i,a in enumerate(seq):
                        pwm[i][a] += seqwt

                for i in range(L):
                    tot = sum( pwm[i].values() )
                    for a in alphabet:
                        pwm[i][a] /= tot
                return pwm


    def _generate_vl_jl(self, chain):
        """
        Parameters
        ----------
        chain : str

        Returns
        -------

        """
        self.vl_nbr, self.jl_nbr = self._get_counts_lists_from_tcr_indices(self.matched_tcrs_plus_nbrs, chain)
        self.vl, self.jl = self._get_counts_lists_from_tcr_indices(self.matched_tcrs, chain)


    def _get_counts_lists_from_tcr_indices(self, indices, chain):
        vcounts = {}
        jcounts = {}
        for ii in indices:
            tcr = self.tcrs[ii]
            if chain in ['A',"G"]:
                vrep,jrep = tcr[6: 8]
            else:
                vrep,jrep = tcr[8:10]
            vcounts[vrep] = vcounts.get(vrep,0)+1
            jcounts[jrep] = jcounts.get(jrep,0)+1
        vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.items()] )
        jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.items()] )
        return self._get_counts_list_condensing_alleles(vstring, self.rep2label_rep, self.rep2label_rep_color),\
               self._get_counts_list_condensing_alleles(jstring, self.rep2label_rep, self.rep2label_rep_color)


    def _get_counts_list_condensing_alleles(self, counts_string, rep2label_rep, rep2label_rep_color ):
        counts ={}
        for tag,count in [x.split(':') for x in counts_string.split(',') ]:
            rc = ( rep2label_rep[ tag ][4:], rep2label_rep_color[ tag ] )
            counts[rc] = counts.get(rc,0)+float(count)
        return [ (y,x[0],x[1]) for x,y in counts.items() ]
