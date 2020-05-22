from . import all_genes_db
from collections import Counter
import numpy as np
from tcrdist.pairwise import hm_metric, nw_metric
import warnings

class TCRcodon():
    """  
    TCRcodon is a class for imputing nucleic acid sequences 
    from cdr3 amino acids and predicted gene usage.
    In particular, TCRcodon helps when only amino acid cdr3 d
    ata and v-gene usage data are available.

    TCRcodon assumes that cdr3 amino acid seqs
    begin with (C) Cystine and ends with (F) Phenylalanine.
    
    Methods
    -------
    guess_reverse_translation : func

    get_best_j_gene : func

    _calculate_gene_specific_codon_usage : func

    Attributes
    ----------
    db_file : str
        'alphabeta_db.tsv' or 'gammadelta_db.tsv' 
    organism : str
        'human' or 'mouse' 
    all_genes : dict
        all_genes is automatically generated at init
        based on db_file and organism arguments. i.e.,
        all_genes_db.all_genes_db(db_file)[self.organism]
    translation_table : dict
        dictionary with uppercase codons as keys
        (e.g., 'ATG') mapping to values uppercase 
        aminos acid codes (e.g., 'M')
    overall_codon_usage : list
        list of lists of all codons used in V,J,D genes
        of speciefied organism
    codon_usage_frequency : collections.Counter
        Counter dict showning how frequentlya particular
        codon was used in the V,J,D genes
        of speciefied organism
    reverse_translation_table : dict 
        dictionary of tuples. Keys are amino acids (e.g., 'G')
        mapping to sorted tuples e.g., 
        'G' : [('GGA', 774), ('GGG', 454), ('GGC', 383), ('GGT', 254)]
        The first tuple shows the most commonly used codon for 
        that particular amino acid.
    cdr3_v_aa_dict = cdr3_v_aa_dict
        dictionary with keys as cdr3 contribtion 
        of each possible V gene e.g., {'CAVR': 'TRAV7-6*02',...}
    cdr3_j_aa_dict = cdr3_j_aa_dict
        dictionary with keys as cdr3 contribtion 
        of each possible J gene e.g., {'DSGYNKLTF': 'TRAJ11*01',..}  
    gene_specific_translation_tables : dict
        dictionary keys are V and J gene names mapping 
        to the specific amino acids and codons they contribute
        to a CDR3 (e.g., 
        'TRBV12-1*01': [('C', 'TGT'), ('A', 'GCC'), ('S', 'AGC'), ('S', 'TCT'),  ('L', 'CTC')]   
    """
    def __init__(self,
                 db_file = "alphabeta_db.tsv",
                 organism = "mouse" ):
        """
        Parameters
        ----------
        db_file : str
            'alphabeta_db.tsv' or 'gammadelta_db.tsv' 
        organism : str
            'human' or 'mouse' 
        
        Example
        -------
        >>> from tcrdist.reverse_translate import TCRcodon
        >>> tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
        """
        self.db_file = db_file
        self.organism = organism
        self.all_genes = all_genes_db.all_genes_db(db_file)[self.organism]
        self.translation_table = { 
                'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
                'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
                'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
                'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
                'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
                'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
                'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
                'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
                'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
                'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
                'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
                'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
                'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
                'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
                'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
                'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
            }
        self.overall_codon_usage = None
        self.codon_usage_frequency = None
        self.cdr3_v_aa = None
        self.cdr3_v_codon_usage = None
        self.cdr3_j_aa = None
        self.cdr3_j_codon_usage = None
        self.reverse_translation_table = None
        self.gene_specific_translation_tables = None

        # THIS IS THE MAJOR INITIALIZATION FUNCTOIN
        self._calculate_gene_specific_codon_usage()
    
    def get_best_j_gene(self, aa_seq, cdr3_ref_dict = None, verbose = False):
        """
        WARNING: THIS SHOULD NOT BE RELIED ON FOR ANYTHING AT THIS POINT 
        BUT IT DOES FILL MISSING J-GENES WITH A 
        GUESS IF ALL YOU HAVE IS AMINO ACID SEQS.

        Parameters
        ----------
        aa_seq : str
            string of uppercase amino acid letters
        cdr3_ref_dict : dict
            defaults to self.cdr3_j_aa_dict
        verbose : bool
            defaults ot False, 

        Returns
        -------
        top_hit : str
            gene_name of the 'best' j-gene hit

        Example 
        -------
        >>> tc.get_best_j_gene(aa_seq =  'CALGELEVVRGYTYTDKLIF', verbose = False)
        'TRDJ1*01'

        Notes
        -----
        TODO: this is an experimental feature using nw_metric, but differences
        in length of CDR3 contribution from different genes may 
        introduce a massive bias. This is better than nothing for 
        making reasonable codon guess when no other information is present
        but does not currently perform very well. A new metric may improve 
        this to a limited extent.    
        """
        warnings.warn("WARNING: get_best_j_gene() IS EXPERIMENTAL AND"\
        "SHOULD NOT BE RELIED ON FOR ANYTHING AT THIS POINT "\
        "BUT IT DOES FILL MISSING J-GENES WITH A DECENT GUESS "\
        "IF ALL YOU HAVE IS AMINO ACID SEQS.")
        
        if cdr3_ref_dict is None:
            cdr3_ref_dict = self.cdr3_j_aa_dict
        
        hits = [(name, nw_metric(aa_seq[-len(k)::], k), k, aa_seq[-len(k)::],aa_seq) for k,name in cdr3_ref_dict.items()]
        top_hit = self._sort_tuple(hits, pos = 1, reverse = False)[0][0]
        if verbose:
            warnings.warn("-------")
            warnings.warn(",".join(map(str,self._sort_tuple(hits, pos = 1, reverse = False))))
            warnings.warn("-------")
            warnings.warn(top_hit)
            warnings.warn("-------")
        
        return top_hit

    
    def guess_reverse_translation(  self,
                                    v_gene_name, 
                                    j_gene_name, 
                                    cdr3_aa,
                                    f = None,
                                    verbose = False):
        """
        Given v and j gene names and cdr_aa amino acid (C....F),
        make an inteligent guess at the actual nucleic acid
        codons used. 

        Parameters
        ----------
        v_gene_name : str
            string name of v gene e.g., 'TRBV29*01' 
        j_gene_name : str
            string name of j gene e.g., 'TRBJ1-5*01' 
        cdr3_aa : str
            string of uppercase amino acid characters 'CASSEGEAPLF'
        f : function
            function that replace amino with codon when
            it can't be imputed from the v or j gene.
            This defaults to TCRcodon.reverse_trans_func
        verbose : bool
            defauls to False. If True, will output diagnostic
            stepwise imputation results starting with (V)
            codons and then (J) codons.
        
        Example
        -------
        >>> from tcrdist.reverse_translate import TCRcodon  
        >>> tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
        >>> tc.guess_reverse_translation(v_gene_name= 'TRBV29*01' , j_gene_name= 'TRBJ2-2*01' , cdr3_aa = 'CASSPTGQLYF')
        'TGTGCTAGCAGTCCTACCGGGCAGCTCTACTTT'
        """
        # <f> is the default function for imputing a codon from a amino acid
        # default uses the most frequent codon for a given amino acid 
        # in all organism specific V,J, and D genes
        if f is None:
            f = self.reverse_trans_func
    
        v_translation = list(self.gene_specific_translation_tables[v_gene_name])
        j_translation = list(self.gene_specific_translation_tables[j_gene_name])
        cdr3_aa_list = [aa for aa in cdr3_aa]
        if verbose: 
            warnings.warn(f"target: {cdr3_aa_list}")
            warnings.warn(f"{v_gene_name}:{v_translation}")
            warnings.warn(f"{j_gene_name}:{j_translation}")
        
        forward_translation = []
        backward_translation = []
        # v_gene_pass 
        
        for amino_acid in cdr3_aa_list:
            if v_translation:
                aa, codon = v_translation.pop(0)
                if amino_acid == aa:
                    forward_translation.append(codon)
                else:
                    forward_translation.append("XXX")
            else:
                forward_translation.append("XXX")
        
        for amino_acid in cdr3_aa_list[::-1]:
            if j_translation:
                aa, codon = j_translation.pop(-1)

                if amino_acid == aa:
                    backward_translation.insert(0,codon)
                else:
                    backward_translation.insert(0,"XXX")
            else:
                backward_translation.insert(0,"XXX")

        combined_translation = [v if v != "XXX" else j for v,j in zip(forward_translation, backward_translation)]
        combined_translation_with_unknowns = [codon if codon != "XXX" else f(aa) for codon, aa in zip(combined_translation, cdr3_aa_list )]   
        if verbose: 
            warnings.warn(f"forward_translation  {forward_translation}")
            warnings.warn(f"backward_translation {backward_translation}")
            warnings.warn(f"combined_translation {combined_translation }")
            warnings.warn(f"with_unknowns        {combined_translation_with_unknowns}")
        
        return "".join(combined_translation_with_unknowns)

    def _sort_tuple(self,tup, pos =1, reverse = True): 
        """
        Parameters
        ----------
        tup : list of tuples
            list of tuples of any fixed length
        pos : int
            index of the tuple position used for sorting
        reverse : bool
            If True largest tuples come first
        
        Returns
        -------
        sorted_tup : list
            list of tuples

        Example
        -------
        >>> from tcrdist.reverse_translate import TCRcodon  
        >>> tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
        >>> tc._sort_tuple([('A',10),('B',20), ('C', 15)])
        [('B', 20), ('C', 15), ('A', 10)]
        """
        sorted_tup = sorted(tup, key = lambda x: x[pos], reverse = reverse)
        return sorted_tup

    def _count3(self,nucseq, protseq, nucseq_offset):
        """
        Parameters
        ----------
        nucseq : str
            nucleic acid string
        protseq : str
            corresponding amino acid string
        nucseq_offset : int
            which nucleic acid in string to start with to match protseq (usually zero)
        
        Returns
        -------
        info : tuple 
            3-part tuple (codons : list, amino_acids : list, translation : list)
        
        Notes
        -----
            if everything is working amino_acids and translation 
            components of results should match in so much as 
            np.all(amino_acids == translation[0:len(amino_acids)])

        Example
        -------
        >>> from tcrdist.reverse_translate import TCRcodon  
        >>> tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
        >>> tc._count3('ATGATG', 'MM', 0)
         (['ATG', 'ATG'], ['M', 'M'], ['M', 'M'])
        """
        nucseq = nucseq[nucseq_offset:]
        codons = [nucseq[i:i+3].upper() for i in range(0, len(nucseq), 3)]
        amino_acids = [aa for aa in protseq]
        translation = list()
        for codon in codons:
            try: 
                translation.append(self.translation_table[codon])
            except KeyError:
                translation.append("X")
        
        info = (codons, amino_acids, translation)
        return info 

    def _calculate_gene_specific_codon_usage(self):
        """
        This is called at initialization of TCRcodon
        """
        # Initialize storage objects
        overall_codon_usage = list()
        cdr3_v_aa = list()
        cdr3_v_aa_dict = dict()
        cdr3_v_codon_usage = list()
        cdr3_j_aa = list()
        cdr3_j_aa_dict = dict()
        cdr3_j_codon_usage = list()
        gene_specific_translation_tables = dict()

        # Note: self.all_genes is organism and alphabeta/gammadelta db specific
        # Note: we are going to loop through all the genes in reference database
        # Note: There are two GOALS here:
        #   1. find the most frquently used codon for a given amino acid in this reference set
        #   2. generate a gene_specific_translation_tables dictionary
        #       that we can use to reverse translate the ends of a given 
        #       cdr3, when we know the germline genes that generated it
        for gene_name, tcr_gene_object in self.all_genes.items():
            # Get codons associated with each amino acid
            (codons, amino_acids, translation) = self._count3(tcr_gene_object.nucseq, 
                                                        tcr_gene_object.protseq,
                                                        tcr_gene_object.nucseq_offset)
            # Check match between amino acid and translation
            assert np.all(amino_acids == translation[0:len(amino_acids)])
            
            # Collect V, D, and J raw codons used for goal 1 (see above)
            overall_codon_usage.append(codons[0:len(amino_acids)])
            
            # Handle V genes for goal 2 (see above)
            if tcr_gene_object.region == "V":
                cdr3_string = tcr_gene_object.cdrs[3].replace(".", "")
                cdr3_start  = tcr_gene_object.protseq.find(cdr3_string)
                cdr3_end    = cdr3_start + len(cdr3_string)
                assert tcr_gene_object.protseq[cdr3_start:cdr3_end] == cdr3_string
                assert "".join([self.translation_table[c] for c in codons[cdr3_start:cdr3_end]]) == cdr3_string
                cdr3_codons = codons[cdr3_start:cdr3_end]
                cdr3_v_codon_usage.append(cdr3_codons)
                cdr3_v_aa.append(cdr3_string)
                cdr3_v_aa_dict[cdr3_string] = gene_name
                # Note: <s> is an amino acid, <c> is a codon
                # for each gene we create a list of ordered tuple mapping amino acid to a specific codon
                gene_specific_translation_tables[gene_name] = [(s,c)  for s,c in zip(cdr3_string, cdr3_codons)]

            # Handle V genes for goal 2 (see above)
            if tcr_gene_object.region == "J":  
                cdr3_string = tcr_gene_object.cdrs[0].replace(".", "")
                cdr3_start  = tcr_gene_object.protseq.find(cdr3_string)
                cdr3_end    = cdr3_start + len(cdr3_string)
                assert tcr_gene_object.protseq[cdr3_start:cdr3_end] == cdr3_string
                assert "".join([self.translation_table[c] for c in codons[cdr3_start:cdr3_end]]) == cdr3_string
                cdr3_codons = codons[cdr3_start:cdr3_end]
                cdr3_j_codon_usage.append(cdr3_codons)
                cdr3_j_aa.append(cdr3_string)
                cdr3_j_aa_dict[cdr3_string] = gene_name
                gene_specific_translation_tables[gene_name] = [(s,c) for s,c in zip(cdr3_string, cdr3_codons)]


        # For Goal 1 (see above) using all codon form V,J,D genes
        # flatten list of list 
        overall_codon_usage_list = [item for sublist in overall_codon_usage for item in sublist]
        # count frequency of each codon
        codon_usage_frequency = Counter(overall_codon_usage_list)
        # make a reverse translation table
        reverse_translation_table = dict()
        for codon, count in codon_usage_frequency.items():
            if len(codon) == 3:
                amino = self.translation_table[codon]
                if amino not in reverse_translation_table.keys():
                    reverse_translation_table[amino] = []
                reverse_translation_table[amino].append((codon, count))

        # TODO: BY SORTING, BY THE MOST FREQUENTLY USED CODON WE CAN USE THAT TO
        # IMPUTE A NUCLEIC ACID CODON THAT ARISE FROM VJ OR VDJ RECOMBINATION. 
        # IN THE FUTURE A TABLE COULD BE MADE FROM ACTUAL CDR3 SEQUENCES WHICH WOULD LIKELY BE MORE ACCURATE.
        # MOREOVER, THIS WOULD BE A CASES FOR MACHINE LEARNING TO GAIN PERFORMANCE, BUT LET'S LEAVE THAT FOR ANOTHER DAY. 
        reverse_translation_table = {k:self._sort_tuple(tup) for k,tup in reverse_translation_table.items()}
        
        self.reverse_translation_table = reverse_translation_table
        self.overall_codon_usage = overall_codon_usage 
        self.codon_usage_frequency =codon_usage_frequency
        self.cdr3_v_aa = cdr3_v_aa 
        self.cdr3_v_codon_usage = cdr3_v_codon_usage
        self.cdr3_j_aa = cdr3_j_aa
        self.cdr3_j_codon_usage = cdr3_j_codon_usage
        self.cdr3_v_aa_dict = cdr3_v_aa_dict
        self.cdr3_j_aa_dict = cdr3_j_aa_dict
        self.gene_specific_translation_tables = gene_specific_translation_tables

    def reverse_trans_func(self, aa):
        """
        Using generated reverse translation table return
        the most frequently used codon for a given amino acid
        
        Parameters
        ----------
        aa : str
            amino acid character uppercase e.g., 'M'
        
        Returns
        -------
        codon : str
            three letter nucleic acid codon sequence, e.g., 'ATG'
        Example
        -------
        >>> from tcrdist.reverse_translate import TCRcodon  
        >>> tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
        >>> tc.reverse_trans_func("M") 
        'ATG'
        >>> tc.reverse_trans_func("G") 
        'GGA'
        """

        codon = self.reverse_translation_table[aa.upper()][0][0]
        return codon

