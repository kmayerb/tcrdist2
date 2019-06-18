# April 25, 2019
import tcrdist as td
import pytest
tempSkip = pytest.mark.skip(reason="Temporarily skipping for efficiency.")

# test that td.processNT returns td.objects.TCRChain instance with either BLAST or Parasail Method
def test_processNT_w_Parasail():
    betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
    betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'
    chainParasail = td.processNT('human', 'B', betaNT, betaQuals, use_parasail= True)
    assert(isinstance(chainParasail, td.objects.TCRChain))

@tempSkip
def test_processNT_w_Blast():
    betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
    betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'
    chainBlast    = td.processNT('human', 'B', betaNT, betaQuals, use_parasail= False)
    assert(isinstance(chainBlast, td.objects.TCRChain))

@tempSkip
def test_processNT_returns_same_result_values():
    betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
    betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'
    chainBlast = td.processNT('human', 'B', betaNT, betaQuals, use_parasail=False)
    chainParasail = td.processNT('human', 'B', betaNT, betaQuals, use_parasail=True)
    ls_b = chainBlast.to_list() # had to create to_list() to get all kwargs in TCRChain
    ls_p = chainParasail.to_list()
    assert(ls_b == ls_p)

@tempSkip
def test_processNT_returns_complete_set_of_result():
    betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
    betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'
    chainBlast = td.processNT('human', 'B', betaNT, betaQuals, use_parasail=False)
    chainParasail = td.processNT('human', 'B', betaNT, betaQuals, use_parasail=True)
    ls_b = chainBlast.to_list() # had to create to_list() to get all kwargs in TCRChain
    ls_p = chainParasail.to_list()
    assert(ls_b == ['jb_gene', 'b_good_hits', 'vb_evalue', 'jb_mm', 'vb_rep', 'vb_countreps', 'cdr3b', 'vb_alignlen', 'cdr3b_quals', 'vb_genes', 'jb_countreps', 'vb_bitscore_gap', 'vb_mismatches', 'cdr3b_nucseq', 'jb_alignlen', 'b_status', 'jb_blast_hits', 'vb_gene', 'jb_evalue', 'jb_reps', 'jb_bitscore_gap', 'jb_rep', 'vb_reps', 'jb_mismatches', 'cdr3b_plus', 'vb_blast_hits', 'jb_genes', 'vb_mm'])
    assert(ls_b == ls_p)

@tempSkip
def test_processNT_returns_same_cdr3b_seq():
    betaNT = 'CGGGGGGGGTACCNTTGNTTAGGTCCTCTACACGGTTAACCTGGTCCCCGAACCGAAGGTCAATAGGGCCTGTATACTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAGTCANN'
    betaQuals = '12.12.12.12.12.22.9.8.6.6.6.8.3.0.3.10.3.0.3.10.10.11.20.25.30.37.37.29.27.14.14.15.27.30.41.47.36.50.50.50.42.42.57.57.43.47.53.47.47.47.47.47.47.50.54.57.57.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.57.57.57.57.59.59.59.57.57.57.57.57.57.57.57.59.57.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.68.59.59.59.59.59.57.57.57.59.57.57.43.37.28.28.21.28.23.37.28.30.15.19.17.15.21.20.25.3.0.0'
    chainBlast = td.processNT('human', 'B', betaNT, betaQuals, use_parasail=False)
    chainParasail = td.processNT('human', 'B', betaNT, betaQuals, use_parasail=True)
    assert(chainBlast.__getattr__("cdr3b") == chainParasail.__getattr__("cdr3b"))
