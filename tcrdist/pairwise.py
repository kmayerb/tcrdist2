import parasail
import multiprocessing
from scipy.spatial import distance
import itertools
import random
import numpy as np

def apply_pw_distance_metric_w_multiprocessing(sequences,
                                               f,
                                               processes = multiprocessing.cpu_count()):
    """
    Function that computes pairwise distance between a list of sequences,
    using python's multiprocessing package to parralelize the distance calulation
    to multiple python interpretors based on the available cpu.

    This is the core function of this module!!!!!

    It is flexible, as f can be any function taking as input a list of indices.
    [(i,j),(i,j), ... (i,j)].

    unique_seqs an orderd list containing all of the sequence to be compared
    are put in the global namespace of each interpretor, and f can be defined
    to calculate the distance metric as appropriate.

    For a first example:

    **_f_pwdist_parallel_using_nw_metric** uses  **nw_metric()**

    _f_pwdist_parallel_using_nw_metric():
        output_tuples = []
        for i,j in indices:
            d = nw_metric(unique_seqs[i], unique_seqs[j])
            output_tuples.append((i,j,d, unique_seqs[i], unique_seqs[j]))
        return(output_tuples)

    For a second example:

    **_f_pwdist_parallel_using_distance_wrapper** calls a distance_wrapper()

    output_tuple = []
    for i,j in indices:
        d = distance_wrapper(unique_seqs[i], unique_seqs[j])
        output_tuple.append((i,j,d, unique_seqs[i], unique_seqs[j]))
    return(output_tuple)


    distance_wrapper():

    def distance_wrapper(a,b):
        return(float(SequencePair(a, b).hamming_distance) )

    This currently wraps scipy.hamming distance after aligning sequences


    Parameters
    ----------
    sequences : list
        list of strings containing amino acid letters
    processes : ind
        number of available cpus defaults to multiprocessing.cpu_count()
    f : function
        function that references unique_seqs and can take a list of (i,j) tuples

    Returns
    -------
    multiprocessed_result : list
        list of lists containing five part tuples (int, int, float, str, str)

    Examples
    --------
     pairwise.apply_pw_distance_metric_w_multiprocessing(sequences = sequence
     s, f = pairwise._f_pwdist_parallel_using_distance_wrapper,
     processes = 2)

     pairwise.apply_pw_distance_metric_w_multiprocessing(sequences = sequences,
     f = pairwise._f_pwdist_parallel_using_nw_metric,
     processes = 2)

    """

    error_message = "in a_apply_pairwise_multiprocessing() you can not \
    specify a higher processes number than available CPUs."

    if processes > multiprocessing.cpu_count():
        raise RuntimeError(error_message)

    # this produces chunked indices a list of lists containing (i,j) tuples
    indices = get_chunked_pwdist_indices(sequences = sequences,
                                         processes=processes)

    # sets the sequences to the namespace of each multiprocerss
    def set_global(x):
        unique_seqs = x
        global unique_seqs

    # creates pool, and runs intializer functoin
    p = multiprocessing.Pool(processes = processes,
                             initializer = set_global,
                             initargs=(sequences,))

    # map function to the parralel processes
    multiprocessed_result  = p.map(f, indices)
    p.close()
    p.join()
    return(multiprocessed_result)

def get_random_amino_acid(n):
    """
    Function that results in random amino acid seq of length n
    for testing and benchmarking
    """
    return(''.join([random.choice('GPAVLIMCFYWHKRQNEDST') for i in xrange(n)]))

def get_k_ramdom_amino_acid_of_length_n(k = 1000, n = 20):
    """
    Function that results in a list of k amino acid seqs of length n
    for testing and benchmarking
    """
    return([get_random_amino_acid(n) for i in xrange(k)])


def _partition(l,n):
    """
    Function that takes a list and maximum number of elements,
    and break the list into sublists

    Parameters
    ----------
    l : list
    n : int

    Returns
    -------
    list of lists

    """
    return([l[i:i + n] for i in xrange(0, len(l), n)])


def get_pwdist_indices(sequences):
    """
    From a list of sequences get lower triangle indices tuples (i,j)
    for pairwise comparisons.

    Parameters
    ----------
    sequences : list of strings
        list of (likely amino acid) strings

    Returns
    -------
    ind_tuples : list
        ist of tuples (i,j)

    """
    ind_tuples = []
    L = len(sequences)
    for i, j in itertools.product(range(L), range(L)):
        if i <= j:
            ind_tuples.append((i,j))
    return(ind_tuples)


def get_chunked_pwdist_indices(sequences, processes):
    """
    Function that takes given a list of sequences, and return a
    chunked set of indices equal to the number of processes.

    Parameters
    ----------
    sequences : list of strings
        list of (likely amino acid) strings

    Returns
    -------
    ind_chunks : list
        list of lists containing tuple indices

    Examples
    -------
    get_chunked_pwdist_indices(sequences = ("A","B","C"), processes = 2)
    >>> (1,2), (1,3)],  [(2,3)] ]
    """

    ind = get_pwdist_indices(sequences)
    chunk_size = (len(ind) / processes) + 1
    ind_chunks = _partition(ind, chunk_size)
    return(ind_chunks)

def _f_pwdist_serial_using_nw_metric(i,j):
    """
    Function for comparing (f_s) serial vs parralelel processing speed.
    This function does distance matrix calculation on at a time.

    ! unique_seqs must be in the namespace when this function is called

    Parameters
    ----------
    i : integer
        integer i index
    j : integer
        integer j index

    Returns
    -------
    five part tuples (int, int, float, str, str)

    """
    d = nw_metric(unique_seqs[i], unique_seqs[j])
    return( (i,j,d, unique_seqs[i], unique_seqs[j]))

def _f_pwdist_parallel_using_nw_metric(indices):
    """
    Function for parralel processing of pairwise distance method using
    hard coded call to the nw_metric

    Parameters
    ----------
    indices: list
        list of tuples (i,j) containing int index refs for pairwise comparison

    Returns
    -------
    output_tuples : list
        list of five part tuples (int, int, float, str, str)

    """
    output_tuples = []
    for i,j in indices:
        d = nw_metric(unique_seqs[i], unique_seqs[j])
        output_tuples.append((i,j,d, unique_seqs[i], unique_seqs[j]))
    return(output_tuples)

def _f_pwdist_parallel_using_distance_wrapper(indices):
    """
    Function for parralel processing of pairwise distance method using
    call to the distacne_wrapper which can be configured

    Parameters
    ----------
    x :

    Returns
    -------
    x :


    """
    output_tuples = []
    for i,j in indices:
        d = distance_wrapper(unique_seqs[i], unique_seqs[j])
        output_tuples.append((i,j,d, unique_seqs[i], unique_seqs[j]))
    return(output_tuples)

def _pack_matrix(chunked_results):
    """
    Function taking a list of tuple generated from a f_pwdist_parallel_()
    function and generaing a np.matrix.

    Parameters
    ----------
    chunked_results : list
        list of lists of five part tuples (int, int, float, str, str)

    Returns
    -------
    pwdist : np.matrix
        matrix with chunked result applied to i,j position


    """
    indices = [item for sublist in chunked_results for item in sublist]
    matrix_dim = max([x[0] for x in indices])+1
    pwdist = np.nan * np.zeros((matrix_dim, matrix_dim))
    for i,j,d,x,y in indices:
        pwdist[j, i] = d
        pwdist[i, j] = d
    return(pwdist)

def flatten(l):
    return([item for sublist in l for item in sublist])

def nw_metric(s1, s2):
    """
    Function applying Parasail's Needleman Wuncsh Algorithm to get a distance
    between any two sequences.

    May or may not produce a true metric. Details in:
    E. Halpering, J. Buhler, R. Karp, R. Krauthgamer, and B. Westover.
    Detecting protein sequence conservation via metric embeddings.
    Bioinformatics, 19 (sup 1) 2003

    Parameters
    ----------
    s1: string
        string containing amino acid letters

    s2: string
        string containing amino acid letters

    Returns
    -------
    D : float
        distance betwee

    """
    xx = parasail.nw_stats(s1, s1, open=3, extend=3, matrix=parasail.blosum62).score
    yy = parasail.nw_stats(s2, s2, open=3, extend=3, matrix=parasail.blosum62).score
    xy = parasail.nw_stats(s1, s2, open=3, extend=3, matrix=parasail.blosum62).score
    D = xx + yy - 2 * xy
    return D






class SequencePair:
    """
    SequencePair

    class that takes a pair of cdr3 or other sequences, aligns them and
    calculated hamming distance (0 - 1.0)

    !!!!!!!!!!WARNING: TODO: GAP AND EXTEND AND OPENNING PENALTIES ARE NOT
    YET VALIDATED COMPARED TO ORGINAL TCRDist

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
    open_penalty : int
        open penalty for parasail.nw
    extend_penalty : int
        extend penalty for parasail.nw
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
        self.open_penalty= 3
        self.extend_penalty = 3
        self.hamming_distance = None
        self.do_align = False
        self.do_distance = False

        if (len(s1) >= 5) and (len(s2) >= 5):
            self.do_align = True
            self.run_align()

        if (self.do_align):
            # check alignment returned a string
            if isinstance(self.a1, str):
                # if so, THEN check that alignment length is equal to or greater
                # (i.e. gaps) than the longest of the inputs
                if ( len(self.a1) >= max(len(self.s1),len(self.s2)) ):
                    self.do_distance= True
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
        r = self.align(self.s1, self.s2, self.open_penalty, self.extend_penalty)
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

    def align(self, x,y, open_penalty, extend_penalty, matrix = parasail.blosum62, ):
        """
        Function that takes two strings and returns a
        parasail result object. Uses:

        parasail.nw_trace() - Needleman-Wunsch Global Alignment
        parasail.blossum62 - blossum62 substitution matrix

        The blosum62 matrix may warrant scrutiny:
        Selecting the Right Similarity-Scoring Matrix
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/

        Parameters
        ----------
        x : string
            Amino acid string passed to arg s1 in parasail.sw_trace
        y : string
            Amino acid string passed to arg s2 in parasail.sw_trace
        open_penalty: int
            gap opening penalty
        extent_penalty: int
            gap extension penalty
        matrix:
            parasail.blossum62 - blossum62 substitution matrix

        Returns
        -------
        r : parasail.bindings_v2.Result
            Parasail result from sw_trace allignmnet
            r.traceback.query,


        Raises
        ------
        """

        r = parasail.nw_trace(s1=x,
                        s2=y,
                        extend = extend_penalty,
                        open = open_penalty,
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

def distance_wrapper(a,b):
    return(float(SequencePair(a, b).hamming_distance) )
# NOTE THAT THE apply_pairwise distance function below can be used with any
# function that accepts two strings
def select_unique_sequences(sequences):
    unique_seqs = list(dict.fromkeys(sequences)) # preserves order in python 3
    return(unique_seqs)

def apply_pairwise_distance(sequences,
                            pairwise_distance_function = distance_wrapper):
    """
    Function apply_pairwise_distance takes a list of sequences,
    takes the unique set, and applies
    a pariwise_distance_func to all unique combinations. A hash is returned.

    To be efficient:
    (1) compare unique unique_sequences.
    (2) Do the pairwise comparision in only one direction A:B = B:A

    Parameters
    ----------
    sequences : list of strings
        list of amino acid strings
    pairwise_distance_function: function
        pairwise distance function that takes strings as inputs

    Returns
    -------
    storage : dictionary
        of form d[uniuqe_string1][unique_string2] = numeric result of the pairwise_distance_function
    """

    # storage is a dictionary of dictionaries
    unique_seqs = sequences
    storage = dict.fromkeys(unique_seqs)
    storage ={k:{} for k in storage.keys()}

    init_index = 0
    final_index = len(unique_seqs)
    rng = range(init_index, final_index)

    for i in rng:
        for j in rng:
            storage[unique_seqs[i]][unique_seqs[j]] = \
            pairwise_distance_function(unique_seqs[i], unique_seqs[j])
        # updates the start index by 1,updates the for loop range to avoid dups
        init_index = init_index + 1
        rng = range(init_index,final_index)
    return(storage)

def unpack_dd_to_kkv(dd):
    """
    Function that unpacks a dictionary of dictionaries (dd) to
    a key1,key1,value (kkv)

    dd = { 'A': {'A': 1, 'B': 2, 'C': 3},
           'B': {'B': 4, 'C': 5},
           'C': {'C': 6}
         }

    to

    kkv = { key1  = ['A','A','A','B','B','C'],
            key2  = ['A','B','C','B','C','C'],
            value = [1,2,3,4,5,6]
           }

    Parameters
    ----------
    dd : dictionary
        dictionary of dictionaries pointing to a pairwise value

    Returns
    -------
    kkv : dictionary
        dictionary of arrays


    Example
    -------
    example of what can be done with the output:

    sequences = ["CAGQASQGNLIF","CAGQASQGNLIF","CAGQASQGNLIFA",
    "CAGQASQGNLIAA","CAGQASQGNLIFAAA","CAGQASQGNLIFAAAAA"]

    d = pairwise.apply_pairwise_distance(sequences)
    pandas.DataFrame(pairwise.unpack_dd_to_kkv(dd = d)).head
    <bound method DataFrame.head
                     key1               key2    values
    0       CAGQASQGNLIAA      CAGQASQGNLIAA  0.000000
    1       CAGQASQGNLIAA    CAGQASQGNLIFAAA  0.133333
    2       CAGQASQGNLIAA  CAGQASQGNLIFAAAAA  0.235294
    3        CAGQASQGNLIF      CAGQASQGNLIAA  0.153846
    5       CAGQASQGNLIAA  CAGQASQGNLIFAAAAA  0.235294
    """
    key1 = []
    key2 = []
    value = []
    for i in sorted(dd):
        for j in sorted(dd[i]):
            key1.append(i)
            key2.append(j)
            value.append(dd[i][j])
    kkv = {"key1":key1,"key2":key2,"value":value}
    return(kkv)




def f(index):#= ["CAGQASQGNLIF","CAGQASQGNLIFA","CAGQASQGNLIAA","CAGQASQGNLIFAAA","CAGQASQGNLIFAAAAA", "CAGQASQGNLIFG"]):
    # index is a tuple, first is the rows to loop over,
    # second is the columns
    index_p1 = index[0]
    index_p2 = index[1]
        # storage is a dictionary of dictionaries
    unique_seqs_p1 =[unique_seqs[i] for i in index_p1]
    storage = dict.fromkeys(unique_seqs_p1)
    storage ={k:{} for k in storage.keys()}
    for i in index_p1:
        for j in index_p2:
            storage[unique_seqs[i]][unique_seqs[j]] = \
            distance_wrapper(unique_seqs[i], unique_seqs[j])
            # updates the start index removing first entry, to avoid duplicate calcs
        if(len(index_p2)>1):
            index_p2.pop(0)
    return(storage)



def apply_pairwise_distance_multiprocessing(sequences,
                                            processes = multiprocessing.cpu_count(),
                                            pairwise_distance_function = distance_wrapper):

    """
    Function apply_pairwise_distance takes a list of sequences,
    takes the unique set, and applies a pariwise_distance_func to all unique combinations.
    A hash is returned.

    This version of the function is explicitly set up to chunk indexing
    to allow efficient multi-processing.

    To be efficient:
    (1) compare unique unique_sequences.
    (2) Do the pairwise comparision in only one direction A:B = B:A

     Chunk indices for passing to parralel processes
    (example for 20x20 comparison in 10 chunks)
    This could be more efficient if the first chunks had fewer indices, but
    on large datasets if chunks is set > processes, this should be fine,
    assuming nodes that finish early can take on a new task.
    ([0, 1], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    ([2, 3], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    ([4, 5], [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    ([6, 7], [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    ([8, 9], [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    ([10, 11], [10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    ([12, 13], [12, 13, 14, 15, 16, 17, 18, 19])
    ([14, 15], [14, 15, 16, 17, 18, 19])
    ([16, 17], [16, 17, 18, 19])
    ([18, 19], [18, 19])

    Parameters
    ----------
    unique_seqs : list of strings
        list of amino acid strings
    processes : int
        process corresponding to parrallel process <= of available cores
    pairwise_distance_function: function
        pairwise distance function that takes strings as inputs

    Returns
    -------
    pooled_storage : list of (dictionary of dictionaries)

    """
    error_message = "in apply_pairwise_distance_mulitprocessing() you can not \
    specify a higher processes number than available CPUs."

    if processes > multiprocessing.cpu_count():
        raise RuntimeError(error_message)

    # unique_seq is a list of unique sequences


    chunks = 2*processes
    total_number_of_seqs = len(sequences)
    index = range(0,total_number_of_seqs)
    partition_size = total_number_of_seqs / chunks
    #assert total_number_of_seqs % chunks == 0, "chunk size not a multiple of unique seqs"
    # index is chunked into two parts (row, and column) so that
    # the pairwise distance task can be mapped to multiple processes
    def partition(l,n):
        return([l[i:i + n] for i in xrange(0, len(l), n)])

    def set_global(x):
        unique_seqs = x
        global unique_seqs

    index_p1 = partition(index, partition_size)
    index_p2 = [range(i[0], total_number_of_seqs) for i in index_p1]
    # part1 and part2 indices are passed as a tuple to each process,
    # memory and unique_sequence are shared in memory
    index = [(index_p1[i],index_p2[i]) for i in range(chunks)]
    # Set up pool multiprocessing

    p = multiprocessing.Pool(processes = processes,
                             initializer = set_global,
                             initargs=(sequences,))
    # map function to the parralel processes
    pooled_storage  = p.map(f, index)


    #pooled_storage = parmap.starmap(f, index = index, unique_seqs = sequences)
    p.close()
    p.join()
    return(pooled_storage)





def unpack_pooled_dd_to_kkv(pooled_dd):
    """

    Function that applies unpack_dd_to_kkv() on a list of
    dictionaries

            pooled_dd = [

            { 'A': {'A': 1, 'B': 2, 'C': 3},
              'B': {'B': 4, 'C': 5},
              'C': {'C': 6}},
            { 'D': {'D': 7, 'E': 8, 'F': 9},
              'E': {'E': 10, 'F': 11},
              'F': {'F': 12}}
            ]

            to

            kkvs = [
            { 'key1'  : ['A','A','A','B','B','C'],
              'key2'  : ['A','B','C','B','C','C'],
              'value' : [1,2,3,4,5,6]},
            { 'key1'  : ['D','D','D','E','E','F'],
              'key2'  : ['D','E','F','E','F','F'],
              'value' : [7,8,9,10,11,12]}
            ]

            finally to

            kkv  = {
            'key1'  : ['A','A','A','B','B','C','D','D','D','E','E','F'],
            'key2'  : ['A','B','C','B','C','C','D','E','F','E','F','F'],
            'value' : [1,2,3,4,5,6, 7,8,9,10,11,12]
            }



    Parameters
    ----------
    pooled_dd : list
        list of dicitionary of dictionarie

    Returns
    ----------
    kkv : dictionary
        dictionary of lists (key1, key2, value format)

    Examples
    ----------
    apply_pairwise_distance_mulitprocessing()

    """
    kkvs = [unpack_dd_to_kkv(pooled_dd[i]) for i in range(0,len(pooled_dd))]
    flatten = lambda l: [item for sublist in l for item in sublist]
    kkv = {}
    for key in kkvs[0]:
        kkv[key] =  flatten([x[key] for x in kkvs])

    #kkv =  {"key1" : flatten([x["key1"] for x in kkvs]),
    #        "key2" : flatten([x["key2"] for x in kkvs]),
    #            "value" : flatten([x["value"] for x in kkvs])}

    return(kkv)