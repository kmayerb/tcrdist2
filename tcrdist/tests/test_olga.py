
"""
This is was used when we first incorporated olga but more important tests are at tests/test_pgen.py
"""

import unittest
from tcrdist import olga
import parasail
import pandas as pd
import numpy as np

import tcrdist.olga.load_model as load_model
import tcrdist.olga.generation_probability as generation_probability
import tcrdist.olga.generation_probability as pgen
import tcrdist.olga.sequence_generation as seq_gen
from tcrdist.olga.utils import nt2aa, determine_seq_type
import os.path as op
from tcrdist.paths import path_to_olga_default_models


class test_olga(unittest.TestCase):







    #def test_olga____is_active():
        #[example_df['cdr3_b_aa']


    #def test_something(self):
        #x = True
        #self.assertTrue(x)





    def test_olga_pgen(self, chain_folder = 'human_T_beta', cdr3 = 'CAWSVAPDRGGYTF', v_b = 'TRBV30*01', v_j ='TRBJ1-2*01'):
        """
        NOT AN ACTUAL UNIT TEST, JUST CODE USED TO TEST IF OLGA WORKED AT ALL IN PYTHON 3
        """

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
        self.assertTrue(np.isclose(x,1.203646865765782e-10))



    def test_olga_sgen(self, chain_folder = 'human_T_beta'):
        """
        NOT AN ACTUAL UNIT TEST, JUST CODE USED TO TEST IF OLGA WORKED AT ALL IN PYTHON 3
        """

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
        self.assertTrue(x)


    def test_olga(self,    chain_folder = 'human_T_beta'):
        """
        NOT AN ACTUAL UNIT TEST, JUST CODE USED TO TEST IF OLGA WORKED AT ALL IN PYTHON 3
        """
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



if __name__ == '__main__':
    unittest.main()
