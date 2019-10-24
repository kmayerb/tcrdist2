class StoreIO():
    """
    StoreIO

    Parent class for passing multiple objects between tcrdist2 functions.
    Attributes are init, set, or modified with keyword arguments.
    """
    def __init__(self, **kwargs):
        self.name = 'StoreIO'
        self.valid_attrs = ['a','b','c','d']
        self.valid_attrs_type = [int,int,int,int]
        # sets all valid attributes to None
        [setattr(self, k, None) for k in self.valid_attrs ]
        # updates attributes based on supplied keyword arguments
        [setattr(self, k, v) for k, v in kwargs.items() if k in self.valid_attrs]

    def __str__(self):
        """
        Enables print(StoreIO) to reveal the state of all
        attributes with the SelfStore instance.
        """
        self_attrs = vars(self)
        self_types = [type(v) for k,v in self_attrs.items() if k in self.valid_attrs]
        header = "{} Attributes:\n\t".format(self.name)
        msg = "\n\t".join(["{} : {} : {} : {}".format(x[0],
                                                      type(x[1]),
                                                      hex(id(x[1])),
                                                      x[1])\
                           for x in self_attrs.items()])
        print(header  + msg)

    def set_attrs_with_kwargs(self, validate_kwargs = True, **kwargs):
        """
        Sets attributes based on keyword arguments

        Parameters
        ----------
        validate_kwargs : bool
            if True, attribute will only be added if its keyword is present in
            self.valid_keys
        """
        for key, value in kwargs.items():
            if validate_kwargs:
                if key in self.valid_attrs:
                    setattr(self, key, value)
            else:
                setattr(self, key, value)

    def _validate_attrs(self):
        """
        Calls self._type_check on all pre-specified valid attributes
        """
        [self._type_check(attr_name, attr_type) for attr_name, attr_type\
         in zip(self.valid_attrs, self.valid_attrs_type)]

        return True

    def _type_check(self, attr_name, attr_correct_type):
        """
        Checks that an attribute of if the correct type

        Parameters
        ----------
        attr_name : str
            attribute name
        attr_correct_type : type
            attribute type

        Returns
        -------
        boolean
            True if attribute type is correct

        Raises
        ------
        TypeError if self.attr is does not match correct type

        """
        attr = getattr(self, attr_name)

        if attr is None:
            return True

        if not isinstance(attr, attr_correct_type):
            raise TypeError("{}.{} must be of type: {}".\
                            format(self.name, attr_name, attr_correct_type))
        else:
            return True

    def _coerce_attributes(self):
        """
        Attempts to coerce all attributes to valid type with calls to
        self._type_coerce

        Returns
        -------
        boolean
            True if all coercions occurred without exception being raised
        """
        for attr_name, attr_correct_type in zip(self.valid_attrs, self.valid_attrs_type):
            if getattr(self, attr_name) is not None:
                self._type_coerce(attr_name, attr_correct_type)
        return True

    def _type_coerce(self, attr_name, attr_correct_type):
        """
        Parameters
        ----------
        attr_name : str
        attr_correct_type : type

        Assigns
        -------
        self.attr_name is set as correct type

        Returns
        -------
        boolean
            True if coercion occurred without exception being raised

        Raises
        ------
        ValueError if trying to coerce float to int
        ValueError if attribute cannot be coerced for any other reason

        Notes
        -----
        This is intended primarily for conversion of strings to int, float, list
        """
        if not hasattr(self, attr_name):
            return False
        else:
            obj = getattr(self, attr_name)

            if isinstance(obj, float) and attr_correct_type is int:
                raise ValueError("Cannot coerce {}.{} to {}. Coercion of float to int is not allowed by {}._type_coerce".\
                                 format(self.name, attr_name, attr_correct_type, self.name))

            try:
                obj_coerced_to_correct_type = attr_correct_type(obj)
                setattr(self, attr_name, obj_coerced_to_correct_type)
                return True
            except ValueError:
                raise ValueError("Cannot coerce {}.{} to {}".\
                                 format(self.name, attr_name, attr_correct_type))


# Create a Storage Class with names matching the dataframe
class StoreIOMotif(StoreIO):
    """
    StoreIOMotif(StoreIO)

    Attributes
    ----------
    file_type      : str
        string indicating this came from a MOTIF file
    count          : int
        ? Phil Bradley Please Describe
    expect_random  : float
        ? Phil Bradley Please Describe
    expect_nextgen : float
        ? Phil Bradley Please Describe
    chi_squared    : float
        statistic with higher values indicating enrichment of motif in TCR
        subset relative to the next-gen non-epitope specific reference
    nfixed         : int
        ? Phil Bradley Please Describe
    showmotif      : str
        e.g, 'ALG.G...kvsf' to be ['A', 'L', 'G', '.', 'G', '.', '.', '.', 'k', 'v', 's', 'f']
    num            : int
        ? Phil Bradley Please Describe - number of motif in motif list
    othernum       : int
        ? Phil Bradley Please Describe
    overlap        : int
        ? Phil Bradley Please Describe - the number of TCRs overlapping with the motif
    ep             : str
        epitope recognized by the subset
    ab             : str
        chain of the TCR for the motif
    nseqs          : int
        number of seqs or clone sequences in the TCR subset
    v_rep_counts   : str
        e.g, 'TRAV6D-6*01:23,TRAV6D-6*03:1'
    j_rep_counts   : str
        e.g., 'TRAJ53*01:21,TRAJ56*01:2,TRAJ9*02:1'


    Note
    ----
    The primary attributes of StoreIOMotif come from a motifs file or DataFrame.
    See tcrdist.TCRMotif.find_cdr3()

    These attributes are added by tcrdist.subset.analyze_motif()

    showmotif              : list
    vl_nbr                 : list
    jl_nbr                 : list
    vl                     : list
    jl                     : list
    matches                : list
    nbr_matches            : list
    matched_tcrs_plus_nbrs : list
    matched_tcrs           : list

    """
    def __init__(self, **kwargs):
        self.name = 'StoreIOMotif'
        v_attrs_motif = [  "file_type",
                           "count",
                           "expect_random",
                           "expect_nextgen",
                           "chi_squared",
                           "nfixed",
                           "showmotif",
                           "num",
                           "othernum",
                           "overlap",
                           "ep",
                           "ab",
                           "nseqs",
                           "v_rep_counts",
                           "j_rep_counts"]
        self.valid_attrs = v_attrs_motif
        self.valid_attrs_type = [str,int,float,float,float,
                                 int,str,int,int,int,
                                 str,str,int,str,str]
        # sets all valid attributes to None
        [setattr(self, k, None) for k in self.valid_attrs ]
        # updates attributes based on supplied keyword arguments
        [setattr(self, k, v) for k, v in kwargs.items() if k in self.valid_attrs]




class StoreIOEntropy(StoreIO):
    """
    StoreIOEntropy(StoreIO)

    Class to recieve the complicated output tuple
    from tcrdist.subset.analyze_matches_using_ngseqs()

    Attributes
    ----------
    ? Phil Bradley Please Describe  

    pwm             : dict
        position-wise matrix
    npwm            : dict
        nucleotide_source for the position-wise-matrix
    ng_lenpwm       : dict
        ng (next-gen) pwm
    ng_fwdpwm       : dict
        ng (next-gen) calculated in the forward direction
    ng_revpwm       : dict
        ng (next-gen) calculated in the reverse direction
    fwdpwm          : dict
        forward position-wise-matrix
    revpwm          : dict
        reverse position-wise-matrix
    scale_by_relent : dict
        scale by relative entropy between rep and next-gen background
    ng_fwdseq_reps  : list
        representatives from next-gen forward
    ng_lenseq_reps  : list
        representatives from next-gen reverse
    num_ng_lenseqs  : int
        number of next-gen seqs
    num_ng_fwdseqs  : int
        number of next-gen forward seqs

    Notes
    -----
    The attributes here are used directly for making entropy-aware motif logos
    """
    def __init__(self, **kwargs):
        self.name = 'StoreIOEntropy'
        v_attrs_entropy = [ "pwm",
                            "npwm",
                            "ng_lenpwm",
                            "ng_fwdpwm",
                            "ng_revpwm",
                            "fwdpwm",
                            "revpwm",
                            "scale_by_relent",
                            "ng_fwdseq_reps",
                            "ng_lenseq_reps",
                            "num_ng_lenseqs",
                            "num_ng_fwdseqs"]
        self.valid_attrs = v_attrs_entropy
        self.valid_attrs_type = [dict, dict, dict, dict, dict, dict, dict, dict,
                                list, list, int, int]
        # sets all valid attributes to None
        [setattr(self, k, None) for k in self.valid_attrs ]
        # updates attributes based on supplied keyword arguments
        [setattr(self, k, v) for k, v in kwargs.items() if k in self.valid_attrs]
