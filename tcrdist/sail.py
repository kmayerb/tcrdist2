from __future__ import division
import numpy as np
import parasail as ps

class ParasailMatch:
    """
    An object that has all the values of BlastMatch but is meant for holding outputs of parasail alignment search
    """
    def __init__(self, query_id, hit_id):
        self.query_id = query_id
        self.hit_id = hit_id
        self.evalue =     "NA" #float( evalue )
        self.identities = "NA" #ident
        self.frame =      'NA'
        self.h_start =    "NA" #hstart ## 0-indexed
        self.h_stop =     "NA" #stop
        self.h_strand =   "NA" #h_strand
        self.h_align =    "NA" #hseq
        self.q_start =    "NA" #qstart
        self.q_stop =     "NA" # qstop
        self.q_strand =   "NA" #q_strand
        self.q_align =    "NA" #
        self.middleseq =  "NA" # middleseq
        self.q2hmap =     "NA" #q2hmap ## 0-indexed numbering wrt to fullseq
        self.valid =      "NA" # True


def dna_reverse_complement(string):
    '''
    provides reverse compliment of a DNA string

    comp_map taken from scikit-bio

    :example:
    >>> dna_reverse_complement("ATG")
    CAT
    '''
    comp_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
        'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N',
        'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'y': 'r', 'r': 'y',
        's': 's', 'w': 'w', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
        'h': 'd', 'v': 'b', 'n': 'n', "-" : "-"}
    complement = "".join(([comp_map[nuc] for nuc in string]))
    return (complement[::-1])


def _q_start(query_seq, q_seq):
    """
    returns the starting pytyon string index of query alignment (q_seq) in the query_sequence (query_seq)

    :param query_seq: string full query sequence
    :param q_seq: string alignment sequence (may contain gaps as "-")
    :return:
    :example:
    >>> _q_start(query_seq = "ATGATG", q_seq = "ATG")
    0
    >>> _q_start(query_seq = "ATGATG", q_seq = "GATG")
    2
    >>> _q_start(query_seq="ATGATG", q_seq="GA-TG")
    2
    """
    q_seq = q_seq.replace("-", "") # remove gaps to get index for original sequence
    q_start = query_seq.find(q_seq)
    return(q_start)

def _q_stop(query_seq, q_seq):
    """
    returns the ending python string index of query alignment (q_seq) in the query sequence (query_seq)

    :param query_seq: string
    :param q_seq: string
    :return: integer
    """
    q_seq = q_seq.replace("-" , "")
    q_start = query_seq.find(q_seq)
    q_stop  = q_start + len(q_seq) - 1
    return(q_stop)

def _h_start(hit_seq, h_seq, h_strand):
    """
    returns the ending python string index of query alignment (h_seq) in the hit reference  sequence (hit_seq)

    :param hit_seq:
    :param h_seq:
    :param h_strand:
    :return: integer
    """
    h_seq = h_seq.replace("_", "")
    h_len = len(h_seq)
    if 1 == h_strand:
        h_start = hit_seq.find(h_seq)
    elif h_strand == -1:
        rc_h_seq = dna_reverse_complement(h_seq)      # when strand -1 , reverse_comp of h_seq gets back to reference direction
        rc_hit_seq = dna_reverse_complement(hit_seq)  # rc -> gets original reference in database
        h_start = rc_hit_seq.find(rc_h_seq)           # find where the alignment starts on reference
        h_start = h_start + h_len - 1                 # but becuase, h_start is at the end of the alignment
    else:
        raise ValueError('h_strand must be 1 or -1')
    return(h_start)

def _h_stop(hit_seq, h_seq, h_strand):
    """
    returns the ending python string index of query alignment (h_seq) in the hit reference  sequence (hit_seq)

    :param hit_seq:
    :param h_seq:
    :param h_strand:
    :return: integer

    """
    h_seq = h_seq.replace("_", "")
    h_len = len(h_seq)
    if h_strand == 1:
        h_start = hit_seq.find(h_seq)
        h_stop = h_start + h_len - 1
    elif h_strand == -1:
        rc_h_seq = dna_reverse_complement(h_seq)
        rc_hit_seq = dna_reverse_complement(hit_seq)
        h_stop = rc_hit_seq.find(rc_h_seq)
    else:
        raise ValueError('h_strand must be 1 or -1')
    return(h_stop)


def _get_ids_by_org_chain_region(organism,
                                 chain,
                                 region,
                                 d):
    '''
    returns ids in dictionary (all_genes) matching chain and region

    :param: organism ('human', 'mouse')
    :param: chain string ('A','B')
    :param: region string ('J', 'V')
    :param: d dictionary containing tcrdist.opjects.TCR_Gene

    :returns: ids list of strings

    :example:
    >>> _get_ids_by_org_chain_region(organism = 'human',
                                 chain = 'A',
                                 region = 'J',
                                 d = all_genes)[0:5]
    ['TRAJ15*01', 'TRAJ15*02', 'TRAJ35*01', 'TRAJ50*01', 'TRAJ24*01']
    '''
    ids = []
    for (id, g) in d[organism].iteritems():
        if g.chain == chain and g.region == region:
            ids.append(id)
    return (ids)


def _get_sequence_tuples_from_ids(ids,
                                  organism,
                                  d):
    '''
    returns list of tuples (id, nucseq, strand) from a  list of ids


    :param: ids list of ids (e.g 'TRAJ15*01')
    :param: organism string ('human', 'mouse')
    :param: d ditionary (all_genes)

    :return: tuple (id, nucseq, strand)
        strand = -1 for reverse complement

    :example:
    >>> _get_sequence_tuples_from_ids(ids = ['TRAJ15*01'])
    [('TRAJ15*01',
  'ccaaccaggcaggaactgctctgatctttgggaagggaaccaccttatcagtgagttcca',
  1),
 ('TRAJ15*01',
  'tggaactcactgataaggtggttcccttcccaaagatcagagcagttcctgcctggttgg',
  -1)]
    '''
    return ([(id, d[organism][id].nucseq, 1) for id in ids] +
            [(id, dna_reverse_complement(d[organism][id].nucseq), -1) for id in ids])




def _create_q2hmap(q_seq, h_seq, q_start, h_start, q_strand=1, h_strand=1):
    '''
    create a query to hit map

    produce the query (q) to hit (h) mapping position in the
    primary query sequence to the position of aligned nucleotide
    in the primary hit reference sequence

    The result is a dictionary mapping the original query sequence to aligned positions in the reference.
    This is made difficult by the potential for gaps in either sides of the pairwise alignment.

    Here is a toy example of an alignment of two sequences with a gap

    POS(10)          1         2
    POS(1) 012345678901234567890123  # this is zero indexed position in the orginal sequence
                >                 <
    QRY:   GGGGGCAAAAAAATTGGAAAAAAA  # This is the actual sequence
    QRY A:      CAAAAAAA-TTGGAAAAAAA # SW Alignment
    HIT A:      CAAAAAAATTTGGAAAAAAA
    HIT:      TTCAAAAAAATTTGGAAAAAAA # This is the hit sequence
    POS(10)             1         2
    POS(1)    0123456789012345678901  # this is zero indexed position in the hit sequence
                >                  <

    :param: q_seq    string of the query hit alignment
    :param: h_seq    string of the hit alignment nucleotide
    :param: q_start  integer start index on the query primary seqeuence
    :param: h_start  integer start index on the hit primary sequence
    :param: q_strand integer (1 or -1)
    :param: h_strand integer (1 or -1)

    :return: q2hmap dictionary

    :example:
    >>> _create_q2hmap(q_seq = 'CAAAAAAA-TTGGAAAAAAA', h_seq = 'CAAAAAAATTTGGAAAAAAA', q_start = 5 , h_start = 2 , q_strand = 1, h_strand = 1)
    {5: (2, 'C'),
     6: (3, 'A'),
     7: (4, 'A'),
     8: (5, 'A'),
     9: (6, 'A'),
     10: (7, 'A'),
     11: (8, 'A'),
     12: (9, 'A'),
     13: (11, 'T'),
     14: (12, 'T'),
     15: (13, 'G'),
     16: (14, 'G'),
     17: (15, 'A'),
     18: (16, 'A'),
     19: (17, 'A'),
     20: (18, 'A'),
     21: (19, 'A'),
     22: (20, 'A'),
     23: (21, 'A')}
    '''
    q2hmap = {}

    i_q = 0  # starting counter for a
    i_h = 0  # starting index b
    gaps = '-'  # gap characters

    for i, q in enumerate(q_seq):
        h = h_seq[i]
        if q not in gaps:
            if h in gaps:
                q2hmap[q_start + i_q] = (-1, '-')  # what to return if there is gap in the hit
            else:
                q2hmap[q_start + i_q] = (h_start + i_h, h)
        if q not in gaps: i_q += q_strand
        if h not in gaps: i_h += h_strand
    return (q2hmap)


def _identities(comp):
    '''
    return percent identity out of 100 based on aligned positions ("|") in the comparision sequence

    :param comp: string containing "|" to designate a identical position
    :return: float

    '''
    if len(comp) > 0:
        ident = round(100*comp.count("|") / len(comp))
    else:
        ident = 0.0
    return(ident)


def _evalue_aproximation(S, lam = 1, k = 1):
    """
    some function for approximating e-values from alignment scores produced by parasail Smith Waterman

    :param parasail_result:
    :return:
    """
    bs = (lam *S - np.log(k)) / np.log(2)
    return(1.0 / ( 2.0 ** bs))

def _all_good_hits_with_scores(hits_scores, max_bit_score_delta_for_good_hits):
    """
    return a list of (id, score, evalue) tuples representing good scores
    based on difference between alignment score for each hit and score of the best hit.

    :param hits_scores: list of tuples [(id,score,evalue),(id,score,evalue)] ranked by score
    :param max_bit_score_delta_for_good_hits: integer
    :return: list of 3-part tuples (str,int,float) for (id,score,evalue)
    """
    top_bit_score = hits_scores[0][1]
    return([x for x in hits_scores if (top_bit_score - x[1]) <= max_bit_score_delta_for_good_hits])

def _evalues_score_gap_tuple(ranked_scores):
    """
    ranked list of 3-part tuples, returns a 2-part tuple with max score and gap to next best socre

    :param ranked_scores:  ranked list of 3-part tuples [(id,score,evalue),...,(id,score,evalue)]
    :return: tuple (float, float) 0: ranked alignment score, 1: the gap between best and next best
    """
    if len(ranked_scores) > 1:  # if there was more than one hit
        # ranked_scores of form [(id, score, evalue),(id, score, evalue)]
        score_gap = ranked_scores[0][1] - ranked_scores[1][1]
        r = (ranked_scores[0][2], score_gap)
    elif len(ranked_scores) == 1:  # for the case where there was only one hit
        r = (ranked_scores[0][2], -99)
    else:
        r = (-99, -99)
    return(r)


def _get_hit_parasail(vj,
                      all_genes,
                      organism,
                      ab,
                      blast_seq,
                      user_matrix = ps.matrix_create("ACGT", 3, -5)):
    """
    The workhorse function that can be inserted into the VJ loop in parse_unpaired_dna_sequence_blastn()
    It returns objects that mimic the returns of functions:
        parse_blast_alignments()
        get_all_hits_with_evalues_and_scores()

    :param vj:
    :param all_genes:
    :param organism:
    :param ab:
    :param blast_seq:
    :param user_matrix:

    :return: 2-part dictionary containing "hits" and "hits_scores
        {"hits" : {"tmp": ParasailMatch.instance}, "hits_scores" : [(id,score,evalue),((id,score,evalue))]}
            "hits" returns a ParasailMatch instance that mimics the attributes of a BlastMatch instance
            "hits_scores" returns a list of 3-part tuples withm, 0: hit_id, 1: alignment score, 2: evalue approximation
    """

    ids = _get_ids_by_org_chain_region(organism=organism,
                                       chain=ab,
                                       region=vj,
                                       d=all_genes)

    # seqs is a list of tuples
    # 0: hit_id,
    # 1: hit_seq (full length reference seq)
    # 2: strand (-1 = rev comp)
    seqs = _get_sequence_tuples_from_ids(ids=ids,
                                         organism=organism,
                                         d=all_genes)
    scores = []
    user_matrix = ps.matrix_create("ACGT", 3, -5)
    for i in range(len(seqs)):
        # smith-waterman alignment implemented using parasail
        s = ps.sw_trace(s1=blast_seq,
                        s2=seqs[i][1],
                        extend=3,
                        open=5,
                        matrix=user_matrix)
        scores.append({'score': s.score,
                       'parasail_result': s,
                       'query_seq': blast_seq,
                       'hit_id': seqs[i][0],
                       'hit_seq': seqs[i][1],
                       'h_strand': seqs[i][2],
                       'q_strand': 1}
                      )

    # sort parasail results from highest to lowest alignment score
    scores = sorted(scores, key=lambda x: x['score'], reverse=True)
    id_score_evalue = [(s['hit_id'],s['score'], _evalue_aproximation(s['score'])) for s in scores]

    # add alignment start positions for the highest scoring alignment
    scores[0]["parasail_result"].get_traceback()
    scores[0]["q_seq"]   = scores[0]["parasail_result"]._traceback.query
    scores[0]["h_seq"]   = scores[0]["parasail_result"]._traceback.ref
    scores[0]["comp"]    = scores[0]["parasail_result"]._traceback.comp
    scores[0]["q_start"] = _q_start(q_seq=scores[0]["q_seq"],
                                    query_seq=blast_seq)

    scores[0]["q_stop"]  = _q_stop(q_seq=scores[0]["q_seq"],
                                   query_seq=blast_seq)

    scores[0]["h_start"] = _h_start(h_seq=scores[0]["h_seq"],
                                    hit_seq=scores[0]["hit_seq"],
                                    h_strand=scores[0]["h_strand"])

    scores[0]["h_stop"]  = _h_stop(h_seq=scores[0]["h_seq"],
                                   hit_seq=scores[0]["hit_seq"],
                                   h_strand=scores[0]["h_strand"])

    scores[0]["identities"] = _identities(scores[0]["comp"])

    # add q2hmap for the highest scoring alignment
    q2hmap = _create_q2hmap(q_seq=scores[0]["q_seq"],
                            h_seq=scores[0]["h_seq"],
                            q_start=scores[0]['q_start'],
                            h_start=scores[0]['h_start'],
                            q_strand=scores[0]['q_strand'],
                            h_strand=scores[0]['h_strand'])

    phony_evalue_must_update_function = _evalue_aproximation(scores[0]['score'])

    # bm2 is going to replace teh BlastMatch instance passed by parse_blast_alignments()
    # we only produce it for the top scoring hit scores[0], but the code is written,
    # so that we could produce ParasailMatches in a loop
    bm2 = ParasailMatch(query_id="tmp",
                        hit_id=scores[0]['hit_id'])
    bm2.evalue = phony_evalue_must_update_function #! SHOULD UPDATE
    bm2.identities = scores[0]["identities"]  # percent identities out of 100
    bm2.h_start    = scores[0]['h_start']     # 0-indexed
    bm2.h_stop     = scores[0]['h_stop']
    bm2.h_strand   = scores[0]['h_strand']
    bm2.h_align    = scores[0]['h_seq']
    bm2.q_start    = scores[0]['h_start']
    bm2.q_stop     = scores[0]['h_stop']
    bm2.q_strand   = scores[0]['q_strand']
    bm2.q_align    = scores[0]['q_seq']
    bm2.middleseq  = scores[0]['comp']
    bm2.q2hmap     = q2hmap                    # q2hmap ## 0-indexed numbering wrt to fullseq
    bm2.valid      = "True"                    # valid IF WHAT?
    bm2.frame      = 'NA'

    # results are meant to mimic the outputs in prior functions from blast version:
        # hits = parse_blast_alignments( blast_tmpfile+'.blast', evalue_threshold, identity_threshold )
        # hits_scores = get_all_hits_with_evalues_and_scores( blast_tmpfile+'.blast' ) ## id,bitscore,evalue
    results = {"hits" : {"tmp": bm2}, "hits_scores" : id_score_evalue}
    return(results)


