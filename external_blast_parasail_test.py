# April 25, 2019
import tcrdist as td
#from tcrdist import sail
#from tcrdist import all_genes
import inspect


def ask_processNT_returns_same_results(organism,
                                       chain,
                                       nt,
                                       quals):
    chainBlast = td.processNT(organism, chain, nt, quals, use_parasail=False)
    chainParasail = td.processNT(organism, chain, nt, quals, use_parasail=True)
    ls_b = chainBlast.to_list()  # had to create to_list() to get all kwargs in TCRChain
    ls_p = chainParasail.to_list()

    for attr in ls_b:
        if attr in ls_b and attr in ls_p:
            test_id = (chainBlast.__getattr__( attr ) == chainParasail.__getattr__( attr))
            if not test_id:
                print("\nFor {} : BLAST != PARASAIL".format(attr))
                print("BLAST:\t{}".format(chainBlast.__getattr__( attr )))
                print("PSAIL:\t{}".format(chainParasail.__getattr__( attr)))
        else:
            print(attr +  "missing from Parasail")






if __name__ == "__main__":

# systematically compare alignments and h2qmaps
    import glob
    import os
    fh = open("testfile.txt", 'r')
    storage = []
    for line in fh:
        line = line.strip().split("\t")
        tmp_storage = {}
        blast_seq = line[1]
        id = line[0]
        sail_result = td.sail._get_hit_parasail(vj="V",
                                                all_genes=td.all_genes,
                                                organism="human",
                                                ab="A",
                                                blast_seq= blast_seq)
        tmp_storage["sail"] = sail_result['hits']['tmp']

        td.blast.parse_unpaired_dna_sequence_blastn(organism = "human",
                                                         ab = "A",
                                                         blast_seq = blast_seq,
                                                         info='',
                                                         try_parasail=False,
                                                         use_parasail=False,
                                                         nocleanup=True,
                                                         hide_nucseq=False,
                                                         extended_cdr3=True,
                                                         return_all_good_hits=True,
                                                         max_bit_score_delta_for_good_hits=50)


        list_of_files = glob.glob('./tmp/*inspectV')  # * means all if need specific format then *.csv
        blast_tmpfile = max(list_of_files, key=os.path.getctime)
        evalue_threshold = 1e-1
        identity_threshold = 20
        hits = td.blast.parse_blast_alignments(blast_tmpfile, evalue_threshold, identity_threshold)
        tmp_storage["blast"] = hits['tmp'][0]

        storage.append(tmp_storage)

    import pickle
    print("PICKLING to storage.p")
    pickle.dump( storage, open( "storage.p", "wb" ) )
    fh.close()








    #
    #
    #
    #
    # a1    = "TACAGCAGGGTTTTGTCTGTGATATACCATCAGAATCCTTACTTTGTGACACATTTGTTTGAGAATCAAAATCGGTGAATAGGCAGACAGACTTGTCACTGGAT" \
    #         "TTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATTTGGTTTAACAGAGAGTTTAGTGCCTTTTCCAAAGATGAGATTTCCTTGGCTTGCTTGCC" \
    #         "CAGCACAGAAGTAGATGCCTACATCACTAGGTATGGATGCTGAGATATTCAGGAAGCTGTCCTTTCTGGTTATACCAAACTGAGCAGTCAGTCTTCCATTTGAGGTCAA" \
    #         "TTCACCAGCCTTATATAAGGCTATCAAGAGGACAGGACCTTCCCCAGGGTCCTGCTTGTACCATAGCCAGGTATTAAATATGCTTGAAGAAGTGCAGTTCATGGAGACATCT" \
    #         "TCTCCTTCCTGGATAAACATAGATTGAGGACTCTGATTCAGCCTGTTGACCACCNAGGNGTGGAGTGATCTCAAGCTGTACACCCTACCTTGGCTCTGGCTCTCTTCTTCTTGGCCG" \
    #         "CCGCCCGGGCTTTCGCTCTGGCNCAGCCTAACTTTTTTAACCAAATGCGAAACCGCCTGCCGGTCCCCCAAG"\
    #
    # a1q  = "6.6.6.6.11.9.9.6.6.6.6.7.7.22.32.37.37.37.37.25.24.18.18.15.17.9.9.9.22.19.34.34.37.37.47.47.47.47.53.53.53.50.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.59.68.68.68.57.57.57.57.57.59.59.68.68.68.68.68.68.57.59.59.59.57.57.57.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.62.62.62.62.62.62.68.68.68.68.68.68.68.68.68.68.62.62.62.62.62.62.62.68.68.68.68.68.68.68.62.62.62.62.62.62.62.68.68.68.68.68.62.62.62.62.62.62.62.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.62.62.62.62.68.62.62.68.68.68.68.68.68.68.62.62.62.62.62.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.59.59.68.68.68.59.59.59.57.57.59.59.59.68.68.68.68.59.59.59.57.68.68.68.68.68.59.59.59.68.68.68.68.59.59.59.59.59.59.59.59.68.68.68.59.59.57.59.57.57.57.38.33.12.12.12.27.36.43.43.24.29.28.17.14.38.38.47.47.57.57.57.57.57.57.57.57.57.59.48.48.48.48.57.57.57.47.57.43.57.47.57.57.57.57.46.46.43.43.43.44.57.48.59.59.59.57.57.57.59.54.54.54.48.43.43.43.43.47.47.54.47.47.43.32.37.31.39.35.47.40.47.48.47.47.47.36.36.36.43.37.24.18.18.12.19.19.33.43.46.39.41.23.23.11.8.3.0.3.8.3.0.3.8.9.10.10.11.10.10.9.9.9.9.10.13.18.17.13.8.8.8.10.8.8.7.8.8.8.8.8.6.6.6.6.7.7.9.10.10.10.10.13.13.13.13.8.8.8.8.8.8.8.8.9.15.15.9.8.8.8.9.10.18.19.25.14.10.8.10.8.8.8.8.8.8.8.8.8.8.8.3.0.3.9.10.12.12.10.8.8.8.10.8.8.11.11.10.9.10.10.10.8.11.14.9.10.8.12.9.11.9.12.10.10.10.10.8.9.8.8.8.7.8.8.9.10.10.8.8.7.7"
    #
    # a1s = td.sail._get_hit_parasail(vj="V",
    #                                  all_genes=td.all_genes,
    #                                  organism="human",
    #                                  ab="A",
    #                                  blast_seq=a1)
    #
    # print(">>>> PARASAIL OUTPUT OF _get_hit_parasail() ")
    # print(a1s['hits']['tmp'].__dict__)
    # print("")
    # print("running  td.processing.parse_unpaired_dna_sequence_blastn(use_parasail=False)")
    # td.blast.parse_unpaired_dna_sequence_blastn(organism = "human",
    #                                                  ab = "A",
    #                                                  blast_seq = a1,
    #                                                  info='',
    #                                                  try_parasail=False,
    #                                                  use_parasail=False,
    #                                                  nocleanup=True,
    #                                                  hide_nucseq=False,
    #                                                  extended_cdr3=True,
    #                                                  return_all_good_hits=True,
    #                                                  max_bit_score_delta_for_good_hits=50)
    #
    #
    # blast_tmpfile  = "./tmp/blasttmpuISdU6.fa.blast.inspectV"
    # evalue_threshold = 1e-1
    # identity_threshold = 20
    # hits = td.blast.parse_blast_alignments(blast_tmpfile, evalue_threshold, identity_threshold)
    # print("")
    # print(">>>> BLAST OUTPUT OF parse_blast_alignments() ")
    # print(hits['tmp'][0].__dict__)
    #







    #
    # # breaks j[_A_frame_mismatch, no_VA_Cpos_blastmatch]
    # a15 = "GTCATACAGATTATGTCTGTGATATACACATCAGAATCCTTACTTTGTGACACATTTGTTTGAGAATCAAAATCGGTGAATAGGCAGACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGGTTCTGGATATTTGGAATGACCGTCAAACTTGTCCCTGTCCCAAAATAGAACTGGTTACCGGTGTCCGCCACAGCACAGAAATAAACGCCTGAGTCTGTGGTCTGGGAAGAGGAAATGTACAATAAGCTGTAGCGTTCCGTAGCGACAGTCGTGGCGCTTAATCTTCCATTCTGTTTTGTCCCTGAGGGAATGTAAAACAGGTTGATGAGCTGTCCCCAAGGGTTTTGATGAAACCACTGCAAATTGTTCACAGAGTCAGAATAAATAGCGCTCTCAGTGAAACTGCATTTGACAACCAAGTTATTGTCTGGAGGACACTACTCAACTTGAATAAACACCTGCTTCCCCCGCCCAGAACGCTCCTGAGATGCCACGCAGGCCACCTTATTGGATGGACCGTCTGATCAAGGTATTCACTAATGTTAACCCATAATACTTATCTTTCTTCATGATTTGCACTCACACTCTGCCTATACAGTAATCCTCCTAT"
    # a15q = "6.6.6.6.6.7.7.6.6.6.6.6.8.11.17.17.32.37.37.21.23.14.16.16.34.34.43.39.34.34.37.37.37.37.40.40.40.47.47.57.57.57.68.68.68.68.68.51.51.50.50.50.50.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.57.59.59.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.57.52.52.52.57.57.57.57.57.52.52.57.59.59.59.57.57.57.57.57.57.57.57.59.59.57.57.57.57.57.68.57.57.57.62.52.52.52.62.62.62.62.62.62.62.62.62.62.62.52.52.62.62.62.62.62.62.62.62.62.35.35.35.44.44.45.45.45.62.62.62.62.62.62.62.52.52.52.62.62.62.52.52.62.62.62.62.62.52.52.52.62.62.62.62.62.62.52.52.52.62.62.62.62.62.62.62.62.62.62.62.62.62.62.62.62.62.62.62.52.52.52.62.62.62.62.52.43.52.52.43.43.43.43.35.43.31.43.43.52.52.62.62.62.43.43.31.31.27.43.43.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.57.52.52.52.57.57.57.57.57.57.57.57.57.57.57.52.52.57.52.43.43.57.57.57.57.59.57.57.57.57.57.52.52.52.52.52.52.52.52.52.52.52.52.57.57.54.54.54.54.52.52.52.52.52.59.57.46.43.43.43.57.57.57.32.21.16.9.9.13.23.32.44.44.48.34.26.15.11.11.24.24.31.27.20.20.28.24.30.12.12.12.19.24.41.47.43.40.47.39.32.31.31.23.21.15.13.16.12.11.11.12.21.23.21.26.13.13.12.17.11.11.11.11.27.32.48.48.43.47.39.39.28.26.12.11.11.13.17.25.11.11.10.18.19.16.15.18.15.13.10.10.10.15.17.15.12.12.12.14.20.30.28.18.15.9.9.9.13.17.13.9.8.8.8.8.8.6.8.17.13.12.10.10.9.10.10.9.8.8.8.8.8.8.8.10.10.9.10.9.9.9.10.15.10.12.10.9.8.8.8.8.8.8.10.9.9.10.10.7.7.8.8.8.10.8.8.6.6.8.8.10.11.11.9.8.8.6.6.8.8.8.10.10.6.6.6.6.8.8.15.8.8.8.8.8.8.8.8.8.8.10.6.6.8.8.8.6.6.6.10.8.8.8.8.8.8.8.6.6.8.7.8.8.8.13.8.8.8.6.6.8.6.6.8.8.9.8.8.9.9.9.9.9.8.8.8.8"
    # #    ask_processNT_returns_same_results("human", "A", a15, a15q)
    # a15s = td.sail._get_hit_parasail(vj="V",
    #                                  all_genes=td.all_genes,
    #                                  organism="human",
    #                                  ab="A",
    #                                  blast_seq=a15)
    #
    #
    #
    # print(">>>> PARASAIL OUTPUT OF _get_hit_parasail() ")
    # print(a15s['hits']['tmp'].__dict__)
    # print("")
    # print("running  td.processing.parse_unpaired_dna_sequence_blastn(use_parasail=False)")
    # td.blast.parse_unpaired_dna_sequence_blastn(organism = "human",
    #                                                  ab = "A",
    #                                                  blast_seq = a15,
    #                                                  info='',
    #                                                  try_parasail=False,
    #                                                  use_parasail=False,
    #                                                  nocleanup=True,
    #                                                  hide_nucseq=False,
    #                                                  extended_cdr3=True,
    #                                                  return_all_good_hits=True,
    #                                                  max_bit_score_delta_for_good_hits=50
    #                                                  )
    #
    #
    # blast_tmpfile  = "./tmp/blasttmpICwlgu.fa.blast.inspectV"
    # evalue_threshold = 1e-1
    # identity_threshold = 20
    # hits = td.blast.parse_blast_alignments(blast_tmpfile, evalue_threshold, identity_threshold)
    # print("")
    # print(">>>> BLAST OUTPUT OF parse_blast_alignments() ")
    # print(hits['tmp'][0].__dict__)
    #









    # an example where there is a tie between
    # vb_rep: BLAST != PARASAIL
    # BLAST:  TRBV12 - 4 * 01
    # PSAIL:  TRBV12 - 3 * 01

#    betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGC' \
#             'TGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
#    betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.' \
#                '36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68'
#    ask_processNT_returns_same_results("human" ,"B", betaNT, betaQuals)

