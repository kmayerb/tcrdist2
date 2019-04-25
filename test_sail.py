if __name__ == "__main__":
    from tcrdist import sail  # sail contains external function for working with parasail
    # test_rev_comp
    assert(sail.dna_reverse_complement("ATG") == "CAT")
    assert(sail.dna_reverse_complement("cat") == "atg")
    assert(sail.dna_reverse_complement("cat-cat") == "atg-atg")
    # test _q_start
    assert(sail._q_start(query_seq = "ATGATG", q_seq = "ATG") == 0)
    assert(sail._q_start(query_seq = "ATGATG", q_seq = "GATG") == 2)
    assert(sail._q_start(query_seq="ATGATG", q_seq="GA-TG") == 2)
    # test _q_end
    assert(sail._q_stop(query_seq="ATGATG", q_seq="ATGA") == 3)
    assert(sail._q_stop(query_seq="ATGATG", q_seq="GA-T") == 4)
    # test _h_start
    assert(sail._h_start(hit_seq="ATGATG", h_seq="ATGA", h_strand= 1) == 0)
    assert(sail._h_start(hit_seq="ATGATG", h_seq="GATG", h_strand= 1) == 2)
    # Position on the original given reverse compliment
    assert(sail._h_start(hit_seq="AAACCG", h_seq="CG", h_strand= -1) == 1)
    assert(sail._h_start(hit_seq="AAACCG", h_seq="CC", h_strand= -1) == 2)
    assert(sail._h_stop(hit_seq="AAACCG", h_seq="CG", h_strand= -1) == 0)
    assert(sail._h_stop(hit_seq="AAACCG", h_seq="CC", h_strand= -1) == 1)
    # test identities
    assert(sail._identities("|||||..|") == 75.0)
    assert(sail._identities("") == 0.0)
    assert(sail._identities("|") == 100.0)
    # test all_good_hits_with_scores
    assert( sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 51, 0.01), ("c", 49, 0.01)], 50) == \
    [('a', 100, 0.01), ('b', 51, 0.01)])
    assert(sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)], 50) == \
            [("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)])
    assert (sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)], 10) == \
            [("a", 100, 0.01)])
    print(sail._all_good_hits_with_scores([("a", 100, 0.01), ("b", 52, 0.01), ("c", 51, 0.01)], 10))

