Mixcr
=====

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Import bulk TCR data previously 
parsed with `mixcr <https://mixcr.readthedocs.io/en/master/>`_ into tcrdist, 
using :py:func:`tcrdist.mixcr.mixcr_to_tcrdist2`. 
We strongly recommend cleaning up this import with the function
:py:func:`tcrdist.mixcr.remove_entries_with_invalid_vgene`. It ensures
that only valid v-gene names are included given your chain of interest.  

Example
-------

.. code-block:: python

    import os
    import numpy as np
    import pandas as pd
    from tcrdist.repertoire import TCRrep
    from tcrdist import mixcr

    clones_fn = os.path.join('tcrdist',
                             'test_files_compact',
                             'SRR5130260.1.test.fastq.output.clns.txt')

    df = mixcr.mixcr_to_tcrdist2(chain = "delta", 
                                 organism = "human", 
                                 clones_fn = clones_fn)
    
    df = mixcr.remove_entries_with_invalid_vgene(df, 
                                                 chain = "delta", 
                                                 organism = "human")
    
    df['subject'] = 'SRR5130260.1'
    
    tr = TCRrep(cell_df = df, 
                organism = "human", 
                chains = ['delta'], 
                db_file='gammadelta_db.tsv')  


.. Tip::

    For additional information on using mixcr, 
    see the section on parsing fastq files 
    at the bottom of this page. 


.. automodule:: tcrdist.mixcr

.. automethod:: tcrdist.mixcr.mixcr_to_tcrdist2

.. automethod:: tcrdist.mixcr.remove_entries_with_invalid_vgene



Mixcr: parsing a fastq
----------------------

This example uses a compact-test data file from gamma/delta T cells.

* SRR5130260.1.test.fastq contains a small sample of data from SRR5130255.1.fastq 

The script below provides an example of running 
`mixcr <https://mixcr.readthedocs.io/en/master/>`_ 
within a docker container  `milaboratory/mixcr:3-imgt 
<https://hub.docker.com/r/milaboratory/mixcr>`_.

.. code-block:: bash

    NAME=SRR5130260.1.test.fastq

    docker run --rm\
        -m 4g \
        -v ~/TCRDIST/tcrdist2/tcrdist/test_files_compact/:/work \
        milaboratory/mixcr:3-imgt \
        -c "mixcr align ${NAME} ${NAME}.vdjca --species hsa; \
        mixcr assemble ${NAME}.vdjca ${NAME}.output.clns; \
        mixcr exportClones ${NAME}.output.clns ${NAME}.output.clns.txt; \
        mixcr exportAlignments ${NAME}.vdjca > ${NAME}.result.txt" 


Running this script produces a set of file:

- tcrdist/test_files_compact/SRR5130260.1.test.fastq.output.clns
- **tcrdist/test_files_compact/SRR5130260.1.test.fastq.output.clns.txt**
- tcrdist/test_files_compact/SRR5130260.1.test.fastq.result.txt
- tcrdist/test_files_compact/SRR5130260.1.test.fastq.vdjca

The example above used the clones file `SRR5130260.1.test.fastq.output.clns.txt`
One can also use  `SRR5130260.1.test.fastq.result.txt`
as input for tcrdist2, using :py:func:`tcrdist.mixcr.mixcr_to_tcrdist2` with the 
seqs_fn argument. If `SRR5130260.1.test.fastq.output.clns.txt` is passed to the 
clones_fn argument of :py:func:`tcrdist.mixcr.mixcr_to_tcrdist2`, 
the DataFrame returned will contain  "clone_id" and "count" columns.
