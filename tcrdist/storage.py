class StoreIO():
    """
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

v_attrs_motif = ["file_type",
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

# Create a Storage Class with names matching the dataframe
class StoreIOMotif(StoreIO):
    """

    """
    def __init__(self, **kwargs):
        self.name = 'StoreIOMotif'
        self.valid_attrs = v_attrs_motif
        # sets all valid attributes to None
        [setattr(self, k, None) for k in self.valid_attrs ]
        # updates attributes based on supplied keyword arguments
        [setattr(self, k, v) for k, v in kwargs.items() if k in self.valid_attrs]


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

class StoreIOEntropy(StoreIO):
    """
    Class to recieve the complicated output tuple
    from tcrdist.subset.analyze_matches_using_ngseqs()

    Notes
    -----
    use positional arguments to initialize

    >>> r = analyze_matches_using_ngseqs()
    >>> for i,t in enumerate([dict, dict, dict, dict, dict, dict, dict, dict, list, list, int, int]):
    ...    assert isinstance(r[i], t)
    >>> es = EntropyStack(**r1)

    """
    def __init__(self, **kwargs):
        self.name = 'StoreIOEntropy'
        self.valid_attrs = v_attrs_entropy
        # sets all valid attributes to None
        [setattr(self, k, None) for k in self.valid_attrs ]
        # updates attributes based on supplied keyword arguments
        [setattr(self, k, v) for k, v in kwargs.items() if k in self.valid_attrs]
