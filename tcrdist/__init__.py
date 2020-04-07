from __future__ import absolute_import
from .version import __version__

import logging
logging.basicConfig(filename='tcrdist.log',
                    level=logging.WARNING,
                    format='[%(asctime)s] [%(levelname)s] [%(name)s: %(lineno)s]\n\t%(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    filemode='w')
logger = logging.getLogger('__init__.py')
logger.debug('Beginning package imports')

#from . import processing
from .hello import *
from .processing import processNT, computeProbs, samplerProb
from .tcr_sampler import alpha_cdr3_protseq_probability, beta_cdr3_protseq_probability
from . import util # changed from utils to util to match tcrdist
from . import plotting
from .objects import TCRClone, TCRChain
from . import datasets
from . import distances
from .all_genes import all_genes
from .sail import *
from .setup_blast import install_blast_to_externals
from . import setup_db
from . import install_test_files
#from .test_resources import adaptive100
from .pairwise import *
from . import mappers
from . import vis_tools
from . import repertoire_db
from . import stats


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
           'vis_tools',
           'stats']
