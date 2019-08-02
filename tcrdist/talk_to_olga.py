import time
import tcrdist.olga.load_model as load_model
import tcrdist.olga.generation_probability as generation_probability
import tcrdist.olga.generation_probability as pgen
import tcrdist.olga.sequence_generation as seq_gen
from tcrdist.olga.utils import nt2aa, determine_seq_type
import os.path as op
from tcrdist.paths import path_to_olga_default_models

import numpy as np


#from optparse import OptionParser

def hello_olga():
    print("OLGA FOLDER ON TRACTOR BEAM")

def test_olga_pgen(chain_folder = 'human_T_beta', cdr3 = 'CAWSVAPDRGGYTF', v_b = 'TRBV30*01', v_j ='TRBJ1-2*01'):
    chain_folder = 'human_T_beta'
    params_file_name = op.join(path_to_olga_default_models,
                               chain_folder,
                               'model_params.txt')
    marginals_file_name = op.join(path_to_olga_default_models,
                                  chain_folder,
                                  'model_marginals.txt')
    V_anchor_pos_file = op.join(path_to_olga_default_models,
                                chain_folder,
                                'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = op.join(path_to_olga_default_models,
                                chain_folder,
                                 'J_gene_CDR3_anchors.csv')

    #Load data
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    #Load model
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file_name)

    #Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
    pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
    #Compute some sequence pgens
    x = pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
    assert(np.isclose(x,1.203646865765782e-10))
    return(x)


def test_olga_sgen(chain_folder = 'human_T_beta'):
    chain_folder = 'human_T_beta'
    params_file_name = op.join(path_to_olga_default_models,
                               chain_folder,
                               'model_params.txt')
    marginals_file_name = op.join(path_to_olga_default_models,
                                  chain_folder,
                                  'model_marginals.txt')
    V_anchor_pos_file = op.join(path_to_olga_default_models,
                                chain_folder,
                                'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = op.join(path_to_olga_default_models,
                                chain_folder,
                                 'J_gene_CDR3_anchors.csv')

    #Load data
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    #Load model
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file_name)

    #Process model/data for sequence generation by instantiating SequenceGenerationVDJ
    seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

    #Generate some random sequences
    x = seq_gen_model.gen_rnd_prod_CDR3()
    return(x)


def test_olga(    chain_folder = 'human_T_beta'):
    chain_folder = 'human_T_beta'
    params_file_name = op.join(path_to_olga_default_models,
                               chain_folder,
                               'model_params.txt')
    marginals_file_name = op.join(path_to_olga_default_models,
                                  chain_folder,
                                  'model_marginals.txt')
    V_anchor_pos_file = op.join(path_to_olga_default_models,
                                chain_folder,
                                'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = op.join(path_to_olga_default_models,
                                chain_folder,
                                 'J_gene_CDR3_anchors.csv')

    #Load data
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
    #Load model
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file_name)

    #Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
    pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)

    #Compute some sequence pgens
    x = pgen_model.compute_regex_CDR3_template_pgen('CASSAX{0,5}SARPEQFF')
    assert(np.isclose(x,6.846877804096558e-10))
    x = pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
    assert(np.isclose(x,1.203646865765782e-10))
    x = pgen_model.compute_nt_CDR3_pgen('TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC')
    assert(np.isclose(x,3.9945642868171824e-14))


    #Process model/data for sequence generation by instantiating SequenceGenerationVDJ
    seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

    #Generate some random sequences
    x = seq_gen_model.gen_rnd_prod_CDR3()
    #assert(x == ('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8))
    assert(isinstance(x[0], str))# == 'TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT' )
    assert(isinstance(x[1], str))# == 'CASSEKRQWESGELFF')

    #('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8)
    x = seq_gen_model.gen_rnd_prod_CDR3()
    assert(isinstance(x[0], str))# == 'TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT' )
    assert(isinstance(x[1], str))# == 'CASSEKRQWESGELFF')

    #assert(x == ('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1))
    #('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1)
    x = seq_gen_model.gen_rnd_prod_CDR3()
    assert(isinstance(x[0], str))# == 'TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT' )
    assert(isinstance(x[1], str))# == 'CASSEKRQWESGELFF')
