#cloneId	cloneCount	cloneFraction	clonalSequence	bestVGene	bestDGene	bestJGene	bestVFamily	bestJFamily	aaSeqCDR3	nSeqCDR3	refPoints

import pytest
import os 
import numpy as np
import pandas as pd
from tcrdist.repertoire import TCRrep
from tcrdist import mixcr

def test_load():
    assert os.path.isfile('tcrdist/test_files/contracting_clones_M_alpha.tsv')
    fn =os.path.join('tcrdist','test_files','contracting_clones_M_alpha.tsv')
    assert os.path.isfile(fn)

def test_columns():
    fn =os.path.join('tcrdist','test_files','contracting_clones_M_alpha.tsv')
    df = pd.read_csv(fn, sep = "\t")
    assert np.all(df.columns == ['Rank', 'Read.count', 'Read.proportion', 'CDR3.nucleotide.sequence',
       'bestVGene', 'bestDGene', 'bestJGene', 'bestVFamily', 'bestJFamily',
       'CDR3.amino.acid.sequence', 'nSeqCDR3', 'refPoints', 'CD',
       'cdr3aa_length'])
    
    minervina_alpha_map = {'Rank' : 'clone_id',
                        'Read.count' : 'count',
                        'Read.proportion' : 'prop',
                        'CDR3.nucleotide.sequence': 'cdr3_a_nucseq',
                        'bestVGene' : 'v_b_gene' , 
                        'bestDGene' : 'd_b_gene' ,
                        'bestJGene' : 'j_b_gene' ,
                        'bestVFamily' : 'bestVFamily',
                        'bestJFamily' :'bestVFamily',
                        'CDR3.amino.acid.sequence' : 'cdr3_a_aa',
                        'nSeqCDR3' : 'n_seq_cdr3',
                        'refPoints' : 'ref_points',
                        'CD' : 'cd',
                        'cdr3aa_length' : 'cdr3_a_aa_len'}

    df = df.rename(columns=minervina_alpha_map)
    assert set(df.columns.to_list()) == set(minervina_alpha_map.values()) 
    dftcr = df[['clone_id','count','v_b_gene','j_b_gene','cdr3_a_aa','cdr3_a_nucseq']]


testfiles = [   ('contracting_clones_M_alpha.tsv', 'alpha'), 
                ('contracting_clones_M_beta.tsv',  'beta'),
                ('contracting_clones_W_alpha.tsv', 'alpha'), 
                ('contracting_clones_W_beta.tsv',  'beta'),
                ('expanding_clones_M_alpha.tsv',  'alpha'),
                ('expanding_clones_M_beta.tsv',  'beta'),
                ('expanding_clones_W_alpha.tsv',  'alpha'),
                ('expanding_clones_W_beta.tsv', 'beta')]
@pytest.mark.parametrize("f, my_chain", testfiles)
def test_convert_minervina_to_mixcr_run_tcrdist(f,my_chain):
    
    fn =os.path.join('tcrdist','test_files', f)
    df = pd.read_csv(fn, sep = "\t")
    df['bestVGene'] = df['bestVGene'].apply(lambda s : s + "*00")
    df['bestJGene'] = df['bestJGene'].apply(lambda s : s + "*00")
    
    map_minervina_to_mixcr = \
        {'Rank':'cloneId',
        'Read.count':'cloneCount',
        'Read.proportion':'cloneFraction',
        'bestVGene': 'allVHitsWithScore',
        'bestDGene': 'allDHitsWithScore',
        'bestJGene':'allJHitsWithScore',
        'CDR3.nucleotide.sequence':'nSeqCDR3',
        'CDR3.amino.acid.sequence':'aaSeqCDR3',
        'refPoints':'refPoints'}

    df = df.rename(columns = map_minervina_to_mixcr)
    # CREATE A FAUX MIXCR OUTPUT
    df.to_csv('dfmix.clns.txt', index = False, sep = "\t")
    # USE TCRDIST2 TOOL FOR PORTING MIXCR OUTPUTS
    dfmix = mixcr.mixcr_to_tcrdist2(chain = my_chain,
                                    organism = "human",
                                    clones_fn = 'dfmix.clns.txt')

    if my_chain == "alpha":
        assert set(dfmix.columns.to_list()) == set(['clone_id', 'count', 'v_a_gene', 'd_a_gene', 'j_a_gene','cdr3_a_nucseq', 'cdr3_a_aa'])
    
    assert df.shape[0] == dfmix.shape[0]
    dfmix = mixcr.remove_entries_with_invalid_vgene(dfmix,
                                                    chain = my_chain,
                                                    organism = "human")
    
    dfmix = mixcr.remove_entries_with_invalid_cdr3(dfmix, chain = my_chain)

    if my_chain == "alpha":
        tr = TCRrep(cell_df = dfmix, organism = "human", chains=['alpha'])
    elif my_chain == "beta":
        tr = TCRrep(cell_df = dfmix, organism = "human", chains=['beta'])

    tr.infer_cdrs_from_v_gene(chain = my_chain, imgt_aligned = True)
    tr.cell_df['subject'] = 'M'
    tr.cell_df['epitope'] = 'X'
    
    if my_chain == "alpha":
        tr.index_cols = ['clone_id', 'subject', 'epitope',
                        'v_a_gene','j_a_gene', 
                        'cdr3_a_aa', 'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                        'cdr3_a_nucseq']
    
    elif my_chain == "beta":
        tr.index_cols = ['clone_id', 'subject', 'epitope',
                        'v_b_gene','j_b_gene', 
                        'cdr3_b_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                        'cdr3_b_nucseq']


    tr.deduplicate()

    if my_chain == "alpha":
        tr._tcrdist_legacy_method_alpha()
        assert isinstance(tr.cdr3_a_aa_pw, np.ndarray)
        assert isinstance(tr.paired_tcrdist, np.ndarray)
    elif my_chain == "beta":
        tr._tcrdist_legacy_method_beta()
        assert isinstance(tr.cdr3_b_aa_pw, np.ndarray)
        assert isinstance(tr.paired_tcrdist, np.ndarray)


def test_combine_betas_and_alphas():
    import pytest
    import os 
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist import mixcr
    import multiprocessing
    

    testfiles = [('contracting_clones_M_alpha.tsv', 'alpha'), 
                ('contracting_clones_M_beta.tsv',  'beta'),
                ('contracting_clones_W_alpha.tsv', 'alpha'), 
                ('contracting_clones_W_beta.tsv',  'beta'),
                ('expanding_clones_M_alpha.tsv',  'alpha'),
                ('expanding_clones_M_beta.tsv',  'beta'),
                ('expanding_clones_W_alpha.tsv',  'alpha'),
                ('expanding_clones_W_beta.tsv', 'beta')]

    betas = []
    alphas = []

    for f,my_chain in testfiles:
        fn =os.path.join('tcrdist','test_files', f)
        df = pd.read_csv(fn, sep = "\t")
        df['bestVGene'] = df['bestVGene'].apply(lambda s : s + "*00")
        df['bestJGene'] = df['bestJGene'].apply(lambda s : s + "*00")

        map_minervina_to_mixcr = \
            {'Rank':'cloneId',
            'Read.count':'cloneCount',
            'Read.proportion':'cloneFraction',
            'bestVGene': 'allVHitsWithScore',
            'bestDGene': 'allDHitsWithScore',
            'bestJGene':'allJHitsWithScore',
            'CDR3.nucleotide.sequence':'nSeqCDR3',
            'CDR3.amino.acid.sequence':'aaSeqCDR3',
            'refPoints':'refPoints'}
        


        df = df.rename(columns = map_minervina_to_mixcr)
        # CREATE A FAUX MIXCR OUTPUT
        df.to_csv('dfmix.clns.txt', index = False, sep = "\t")
        dfmix = mixcr.mixcr_to_tcrdist2(chain = my_chain,
                                organism = "human",
                                clones_fn = 'dfmix.clns.txt')
        dfmix['CD'] = df['CD'].copy()
        dfmix = mixcr.remove_entries_with_invalid_vgene(dfmix,
                                                chain = my_chain,
                                                organism = "human")

        dfmix = mixcr.remove_entries_with_invalid_cdr3(dfmix, chain = my_chain)
        dfmix['source'] = f
        
        if my_chain == "alpha":
            alphas.append(dfmix)
        elif my_chain == "beta":
            betas.append(dfmix)
            
    betas_joined = pd.concat(betas)
    alpha_joined = pd.concat(alphas)

    testsets = [(alpha_joined, 'alpha'), 
                (betas_joined, 'beta')]
    
    tcr_rep_results = dict()
    for dfmix, my_chain in testsets:
        if my_chain == "alpha":
            tr = TCRrep(cell_df = dfmix, organism = "human", chains=['alpha'])
        elif my_chain == "beta":
            tr = TCRrep(cell_df = dfmix, organism = "human", chains=['beta'])

        tr.infer_cdrs_from_v_gene(chain = my_chain, imgt_aligned = True)
        tr.cell_df['subject'] = 'M'
        tr.cell_df['epitope'] = 'X'
        
        if my_chain == "alpha":
            tr.index_cols = ['clone_id', 'subject', 'epitope',
                            'v_a_gene','j_a_gene', 
                            'cdr3_a_aa', 'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',
                            'cdr3_a_nucseq', 'CD','source']
        
        elif my_chain == "beta":
            tr.index_cols = ['clone_id', 'subject', 'epitope',
                            'v_b_gene','j_b_gene', 
                            'cdr3_b_aa', 'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa',
                            'cdr3_b_nucseq', 'CD','source']

        tr.deduplicate()
        my_processes = 10
        if my_processes > multiprocessing.cpu_count():
             my_processes = multiprocessing.cpu_count()

        if my_chain == "alpha":
            tr._tcrdist_legacy_method_alpha(processes = my_processes)
            assert isinstance(tr.cdr3_a_aa_pw, np.ndarray)
            assert isinstance(tr.paired_tcrdist, np.ndarray)
        elif my_chain == "beta":
            tr._tcrdist_legacy_method_beta(processes = my_processes)
            assert isinstance(tr.cdr3_b_aa_pw, np.ndarray)
            assert isinstance(tr.paired_tcrdist, np.ndarray)
        
        tcr_rep_results[my_chain] = tr
        # import pickle
        # with open("tcr_rep_results_for_tonight.p", 'wb') as p:
        #     pickle.dump(tcr_rep_results, p)



