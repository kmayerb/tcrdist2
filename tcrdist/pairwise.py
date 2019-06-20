import parasail
from scipy.spatial import distance

class sequence_pair:
    """
    seq_pair

    class that takes a pair of cdr3A sequences, aligns them and calculated hamming distance (0 - 1.0)


    !!!!!!!!!!WARNING: TODO: HAMMING DISTANCE IS ON ALIGNED SEQUENCES USING
    SMITH - WATERMAN, SO IT MAY IGNORE GAPS AT SEQUENCE ENDS
    AND THE EXTEND AND OPENNING PENALTIES ARE NOT YET VALIDATED COMPARED TO ORGINAL TCRDist

    ...

    Attributes
    ----------
    s1 : string
        raw amino acid string
    s2 : string
        raw amino acid string
    a1 : string
        aligned amino acid string
    a2 : string
        aligned amino acid string
    hamming_distance : double
        hamming distance between aligned strings

    Methods
    -------
    run_align(self)
        calls align, assigns a1 and a2
    run_hamming_distance_on_aligned_strings(self)
        calls hamming_distance_on_aligned_strings assigns aligned_hamming_distance
    align(s1,s2)
        smith-waterman alignment depends on parasail
    hamming_distance_on_aligned_strings(x,y)
        haming distance depends on scipy

    """
    def __init__(self,s1,s2):
        if not isinstance(s1, str) and isinstance(s2, str):
            raise TypeError("s1 and s2 must be strings")

        self.s1 = s1
        self.s2 = s2
        self.a1 = None
        self.a2 = None
        self.hamming_distance = None
        self.do_align = False
        self.do_distance = False

        if (len(s1) >= 5) and (len(s2) >= 5):
            self.do_align = True
            self.run_align()

        if (self.do_align):
            # check alignment returned a string
            if isinstance(self.a1, str):
                # if so, THEN check that alignment length is equal to or greater (i.e. gaps) than the shorter of the two sequences
                if ( len(self.a1) >= min(len(self.s1),len(self.s2)) ):
                    self.do_align = True
                    self.run_hamming_distance_on_aligned_strings()
                # Otherwise Hamming Distance Remains None




    def run_align(self):
        """
        Function calls self.align()

        Assigns
        -------
        Assigns aligned self.s1 (query) to attribute self.a1
        Assigns aligned self.s1 (query) to attribute self.a1
        """
        r = self.align(self.s1,self.s2)
        self.a1 = r.traceback.query
        self.a2 = r.traceback.ref

    def run_hamming_distance_on_aligned_strings(self):
        """
        Function calls self.hamming_distance_on_aligned_strings()

        Assigns
        -------
        Assigns value to self.hamming_distance
        """
        d = self.hamming_distance_on_aligned_strings(x = self.a1, y = self.a2)
        self.hamming_distance = d

    def align(self,x,y):
        """
        Function that takes two strings and returns a
        parasail result object. Currently uses blosum62 default matrix.

        This default selection may warrant scrutiny:
        Selecting the Right Similarity-Scoring Matrix
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/

        Parameters
        ----------
        x : string
            Amino acid string passed to arg s1 in parasail.sw_trace
        y : string
            Amino acid string passed to arg s2 in parasail.sw_trace

        Returns
        -------
        r : parasail.bindings_v2.Result
            Parasail result from sw_trace allignmnet
            r.traceback.query,


        Raises
        ------
        """
        matrix = parasail.blosum62
        r = parasail.sw_trace(s1=x,
                        s2=y,
                        extend = 2,
                        open = 1,
                        matrix = matrix)
        return(r)

    def hamming_distance_on_aligned_strings(self, x, y):
        """
        Function that takes two aligned strings, converts to lists, and
        calculates hamming distance between the two arrays.

        Parameters
        ----------
        x : string
            Aligned amino acid string.
        y : string
            Aligned amino acid string.

        Returns
        -------
        d : double
            The Hamming distance between vectors x and y, between 0 and 1

        Raises
        ------
        AssertionError if len(x) != len(y)
        """
        if (len(x) != len(y)):
            raise AssertionError("string length of {} and {} should be equal \
            if they were correctly aligned".format(x,y))

        d = distance.hamming(u = list(x), v = list(y))
        return(d)
