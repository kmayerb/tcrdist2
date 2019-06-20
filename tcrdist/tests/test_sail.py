from tcrdist import sail  # sail contains external function for working with parasail
from tcrdist import pairwise
'''
test_sail contains basic unit tests for the internal functions defined in sail.py

sail.py was added to enable parasail alignments to replace tcrdists original dependency on external blast calls

'''
# test_rev_comp
def test_dna_reverse_complement():
    assert(sail.dna_reverse_complement("ATG") == "CAT")
def test_dna_reverse_complement2():
    assert(sail.dna_reverse_complement("cat") == "atg")
def test_dna_reverse_complement3():
    assert(sail.dna_reverse_complement("cat-cat") == "atg-atg")

def test_q_start():
    assert(sail._q_start(query_seq = "ATGATG", q_seq = "ATG") == 0)
def test_q_start2():
    assert(sail._q_start(query_seq = "ATGATG", q_seq = "GATG") == 2)
def test_q_start3():
    assert(sail._q_start(query_seq="ATGATG", q_seq="GA-TG") == 2)

# test _q_stop
def test_q_stop():
    assert(sail._q_stop(query_seq="ATGATG", q_seq="ATGA") == 3)
def test_q_stop2():
    assert(sail._q_stop(query_seq="ATGATG", q_seq="GA-T") == 4)

# test _h_start
def test_h_start():
    assert(sail._h_start(hit_seq="ATGATG", h_seq="ATGA", h_strand= 1) == 0)
def test_h_start2():
    assert(sail._h_start(hit_seq="ATGATG", h_seq="GATG", h_strand= 1) == 2)
# Position on the original given reverse compliment
def test_h_start3():
    assert(sail._h_start(hit_seq="AAACCG", h_seq="CG", h_strand= -1) == 1)
def test_h_start4():
    assert(sail._h_start(hit_seq="AAACCG", h_seq="CC", h_strand= -1) == 2)

# MAKE SURE GAPS ARE REMOVED
def test_h_start_with_gap():
    assert(sail._h_start(hit_seq="ATGATG", h_seq="A-TGA", h_strand=1) == 0)
def test_h_start_with_gap2():
    assert(sail._h_start(hit_seq="AAACCG", h_seq="C-C", h_strand= -1) == 2)

def test_h_stop():
    assert(sail._h_stop(hit_seq="AAACCG", h_seq="CG", h_strand= -1) == 0)
def test_h_stop2():
    assert(sail._h_stop(hit_seq="AAACCG", h_seq="CC", h_strand= -1) == 1)

# test identities
def test_idenities():
    assert(sail._identities("|||||..|") == 75.0)
def test_idenities2():
    assert(sail._identities("") == 0.0)
def test_idenities3():
    assert(sail._identities("|") == 100.0)

# test all_good_hits_with_scores
def test_all_good_hits_with_scores():
    assert( sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 51, 0.01), ("c", 49, 0.01)], 50) == \
    [('a', 100, 0.01), ('b', 51, 0.01)])
def test_all_good_hits_with_scores2():
    assert(sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)], 50) == \
    [("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)])
def test_all_good_hits_with_scores3():
    assert (sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)], 10) == \
    [("a", 100, 0.01)])
