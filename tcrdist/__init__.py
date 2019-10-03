
from tcrdist.version import __version__

import logging
logging.basicConfig(filename='tcrdist.log',
                    level=logging.WARNING,
                    format='[%(asctime)s] [%(levelname)s] [%(name)s: %(lineno)s]\n\t%(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    filemode='w')
logger = logging.getLogger('__init__.py')
logger.debug('Beginning package imports')

from tcrdist.hello import *
from tcrdist.processing import processNT, computeProbs, samplerProb
from tcrdist.tcr_sampler import alpha_cdr3_protseq_probability, beta_cdr3_protseq_probability
from tcrdist import processing
from tcrdist import util # changed from utils to util to match tcrdist
from tcrdist import plotting
from tcrdist.objects import TCRClone, TCRChain
from tcrdist import datasets
from tcrdist import distances
from tcrdist.all_genes import all_genes
from tcrdist.sail import *
from tcrdist.setup_blast import install_blast_to_externals
#from .test_resources import adaptive100
from tcrdist.pairwise import *
from tcrdist import mappers
from tcrdist import vis_tools
from tcrdist import repertoire_db
from tcrdist import rep_diff as stats



# from . import embedding  (ImportError: libgfortran.so.1 on linux but not on windows, environmental diff?)

__all__ = ['processing',
           'utils',
           'plotting',
           'datasets',
           'distances',
           'embedding',
           'pairwise',
           'repertoire',
           'repertoire_db',
           'mappers',
           'vis_tools']
