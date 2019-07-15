
from .version import __version__

import logging
logging.basicConfig(filename='tcrdist.log',
                    level=logging.WARNING,
                    format='[%(asctime)s] [%(levelname)s] [%(name)s: %(lineno)s]\n\t%(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    filemode='w')
logger = logging.getLogger('__init__.py')
logger.debug('Beginning package imports')

from .hello import *
from .processing import processNT, computeProbs, samplerProb
from tcr_sampler import alpha_cdr3_protseq_probability, beta_cdr3_protseq_probability
from . import processing
from . import util # changed from utils to util to match tcrdist
from . import plotting
from .objects import TCRClone, TCRChain
from . import datasets
from . import distances
from .all_genes import all_genes
from .sail import *

# from . import embedding  (ImportError: libgfortran.so.1 on linux but not on windows, environmental diff?)

__all__ = ['processing',
           'utils',
           'plotting',
           'datasets',
           'distances',
           'embedding']