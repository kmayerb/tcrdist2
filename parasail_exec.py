import parasail as ps
import tcrdist as td
from tcrdist import sail # sail contains external function for working with parasail


def _get_hit_parasail(vj,
                      all_genes,
                      organism,
                      ab,
                      blast_seq,
                      user_matrix = ps.matrix_create("ACGT", 3, -5)):

    ids = sail._get_ids_by_org_chain_region(organism=organism,
                                       chain=ab,
                                       region=vj,
                                       d=all_genes)

    # seqs is a list of tuples
    # 0: hit_id,
    # 1: hit_seq (full length reference seq)
    # 2: strand (-1 = rev comp)
    seqs = sail._get_sequence_tuples_from_ids(ids=ids,
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
    id_score_evalue = [(s['hit_id'],s['score'], sail._evalue_aproximation(s['score'])) for s in scores]

    # add alignment start positions for the highest scoring alignment
    scores[0]["parasail_result"].get_traceback()
    scores[0]["q_seq"] = scores[0]["parasail_result"]._traceback.query
    scores[0]["h_seq"] = scores[0]["parasail_result"]._traceback.ref
    scores[0]["comp"] = scores[0]["parasail_result"]._traceback.comp
    scores[0]["q_start"] = sail._q_start(q_seq=scores[0]["q_seq"],
                                    query_seq=blast_seq)

    scores[0]["q_stop"] = sail._q_stop(q_seq=scores[0]["q_seq"],
                                  query_seq=blast_seq)

    scores[0]["h_start"] = sail._h_start(h_seq=scores[0]["h_seq"],
                                    hit_seq=scores[0]["hit_seq"],
                                    h_strand=scores[0]["h_strand"])

    scores[0]["h_stop"] = sail._h_stop(h_seq=scores[0]["h_seq"],
                                  hit_seq=scores[0]["hit_seq"],
                                  h_strand=scores[0]["h_strand"])

    scores[0]["identities"] = sail._identities(scores[0]["comp"])

    # add q2hmap for the highest scoring alignment
    q2hmap = sail._create_q2hmap(q_seq=scores[0]["q_seq"],
                            h_seq=scores[0]["h_seq"],
                            q_start=scores[0]['q_start'],
                            h_start=scores[0]['h_start'],
                            q_strand=scores[0]['q_strand'],
                            h_strand=scores[0]['h_strand'])

    phony_evalue_must_update_function = sail._evalue_aproximation(scores[0]['score'])

    # bm2 is going to replace teh BlastMatch instance passed by parse_blast_alignments()
    bm2 = sail.ParasailMatch(query_id="tmp",
                        hit_id=scores[0]['hit_id'])
    bm2.evalue = phony_evalue_must_update_function #!!!!!!!!! MUST UPDATE
    bm2.identities = scores[0]["identities"] # percent identities
    bm2.frame = 'NA'
    bm2.h_start = scores[0]['h_start']  # hstart ## 0-indexed
    bm2.h_stop = scores[0]['h_stop']  # stop
    bm2.h_strand = scores[0]['h_strand']  # h_strand
    bm2.h_align = scores[0]['h_seq']
    bm2.q_start = scores[0]['h_start']
    bm2.q_stop = scores[0]['h_stop']
    bm2.q_strand = scores[0]['h_stop']
    bm2.q_align = scores[0]['q_seq']
    bm2.middleseq = scores[0]['comp']
    bm2.q2hmap = q2hmap  # q2hmap ## 0-indexed numbering wrt to fullseq
    bm2.valid = "True"  # IF WHAT?

    # reults are meant to mimic
    # hits = parse_blast_alignments( blast_tmpfile+'.blast', evalue_threshold, identity_threshold )
    # hits_scores = get_all_hits_with_evalues_and_scores( blast_tmpfile+'.blast' ) ## id,bitscore,evalue
    results = {"hits" : {"tmp": bm2}, "hits_scores" : id_score_evalue}

    return(results)


if __name__ == "__main__":
    all_genes = td.all_genes
    organism = "human"
    ab = "A"
    blast_seq = "TACTGCAAGATTTTGTCTGTGATATACACATCAGAATCCTTACTTTGTGACACATTTGTTTGAGAATCAAAATCGGTGAATAGGCAGACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATTTGGTTTAACAGAGAGTTTAGTGCCTTTTCCAAAGATGAGATTTCCTTGGCTTGCTTGCCCAGCACAGAAGTAGATGCCTACATCACTAGGTATGGATGCTGAGATATTCAGGAAGCTGTCCTTTCTGGTTATACCAAACTGAGCAGTCAGTCTTCCATTTGAGGTCAATTCACCAGCCTTATATAAGGCTATCAAGAGGACAGGACCTTCCCCAGGGTCCTGCTTGTACCATAGCCAGGTATTAAATATGCTTGAAGAAGTGCAGTTCATGGAGACATCTTCTCCTTCCTGGATAAACATAGATCGAGGACTCTGATTCAGCTGTTGACCACCACGG"
    user_matrix = ps.matrix_create("ACGT", 3, -5)
    vj = "V"


    replacement = _get_hit_parasail(vj,
                          all_genes,
                          organism,
                          ab,
                          blast_seq,
                          user_matrix=ps.matrix_create("ACGT", 3, -5)):

    hits = replacement["hits"]
    hits_scores = replacement["hits_scores"]
    print(hits)
    print(hits_scores)
