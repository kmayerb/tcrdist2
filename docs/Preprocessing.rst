Preprocessing
=============

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Inputs
######

The most common T cell receptor data type we anticipate supporting in
in tcrdist2:

* tcrdist - current users of TCRdist may wish to use a clones file.
* VDJDB database - T cell receptors with predicted epitope specificity
* 10X Genomics - paired-chain T cell receptor data
* Adaptive Biotechnologies - single chain bulk receptor sequencing data
* iReceptor - compilations of large T cell receptor sequencing efforts

If you are generating custom T cell receptor sequencing data and need help
using tcrdist2, please get in touch (kmayerbl@fredhutch.org).



Bulk Data Processed with MiXCR
##############################

If you have bulk sequencing data, MiXCR is an excellent option for rapidly
annotating the CDRs and identifying gene V and J gene usage. Great documentation
for MiXCR can be found at `mixcr.readthedocs.io/ <https://mixcr.readthedocs.io/en/master/>`_.

.. code-block::

  mixcr align test_R1.fastq test_R2.fastq test2.vdjca --species hsa
  mixcr assemble test2.vdjca output.clns
  mixcr exportClones output.clns output.clns.txt
  mixcr exportAlignments test2.vdjca test2.vdjca.txt
  more test2.vdjca.txt | cut -d'\t' -f3,5,21,30 > ${res}.v_j_cdr3_n_cdr3_aa.txt


10X Genomics Data Files
#######################




Adaptive Data File
##################

See the 3.4.13 SAMPLE EXPORT section of immunoSEQ `Analyzer Manual <https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf>`_

The standard output is expected to have column headers as shown in the table below.

+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| Reference Manual       | Definition                                                                                                  | Example1              | Example2              |
+========================+=============================================================================================================+=======================+=======================+
| nucleotide             | Full length nucleotide sequence                                                                             | GAGTCGCCCAGA          | ATCCAGCCCTCAGA        |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| aminoAcid              | CDR3 amino acid sequence                                                                                    | CASSLSWYNTGELFF       | CASSLGQTNTEAFF        |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| count                  | The count of sequence reads for the given nucleotide sequence                                               | 13                    | 2                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| frequencyCount         | The percentage of a particular nucleotide sequence in the total sample Total length of the CDR3 region      | 0.003421953           | 5.26E-04              |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| cdr3Length             | The highest degree of specificity obtained when identifying the V gene Name of the V family (if identified) | 45                    | 42                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vMaxResolved           | Name of the V gene if identified, otherwise ‘unresolved’                                                    | TCRBV27-01*01         | TCRBV12               |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vFamilyName            | Number of the V gene allele if identified, otherwise ‘1’                                                    | TCRBV27               | TCRBV12               |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vGeneName              | Name of the potential V families for a sequence if V family identity is unresolved                          | TCRBV27-01            |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vGeneAllele            | Name of the potential V genes for a sequence if V gene identity is unresolved                               | 1                     |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vFamilyTies            | Name of the potential V gene alleles for a sequence if V allele identity is unresolved                      |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vGeneNameTies          | The highest degree of specificity obtained when identifying the D gene Name of the D family (if identified) |                       | TCRBV12-03,TCRBV12-04 |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vGeneAlleleTies        | Name of the D gene if identified, otherwise ‘unresolved’                                                    |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dMaxResolved           | Number of the D gene allele if identified, otherwise ‘1’ or ‘unresolved’ if D gene is unresolved            |                       | TCRBD01-01*01         |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dFamilyName            | Name of the potential D families for a sequence if D family identity is unresolved                          |                       | TCRBD01               |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dGeneName              | Name of the potential D genes for a sequence if D gene identity is unresolved                               |                       | TCRBD01-01            |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dGeneAllele            | Name of the potential D gene alleles for a sequence if D allele identity is unresolved                      |                       | 1                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dFamilyTies            | The highest degree of specificity obtained when identifying the J gene Name of the J family (if identified) | TCRBD01,TCRBD02       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dGeneNameTies          | Name of the J gene if identified, otherwise ‘unresolved’                                                    | TCRBD01-01,TCRBD02-01 |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dGeneAlleleTies        | Number of the J gene allele if identified, otherwise ‘1’                                                    |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jMaxResolved           | Name of the potential J families for a sequence if J family identity is unresolved                          | TCRBJ02-02*01         | TCRBJ01-01*01         |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jFamilyName            | Name of the potential J genes for a sequence if J gene identity is unresolved                               | TCRBJ02               | TCRBJ01               |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jGeneName              | Name of the J gene if identified, otherwise ‘unresolved’                                                    | TCRBJ02-02            | TCRBJ01-01            |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jGeneAllele            | Number of the J gene allele if identified, otherwise ‘1’                                                    | 1                     | 1                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jFamilyTies            | Name of the potential J families for a sequence if J family identity is unresolved                          |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jGeneNameTies          | Name of the potential J genes for a sequence if J gene identity is unresolved                               |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jGeneAlleleTies        | Name of the potential J gene alleles for a sequence if J allele identity is unresolved                      |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vDeletion              | The number of deletions from the 3’ end of the V segment                                                    | 0                     | 6                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| n1Insertion            | The number of insertions in the N1 region                                                                   | 2                     | 4                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| d5Deletion             | The number of deletions from the 5’ end of the D segment                                                    | 0                     | 1                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| d3Deletion             | The number of deletions from the 3’ end of the D segment                                                    | 10                    | 5                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| n2Insertion            | The number of insertions in the N2 region                                                                   | 3                     | 2                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jDeletion              | The number of deletions from the 5’ end of the J segment                                                    | 2                     | 1                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vIndex                 | Distance from the start of the sequence (designated 0) to the conserved cysteine residue in the V motif     | 36                    | 39                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| n1Index                | Distance from 0 to the start of the N1                                                                      | 53                    | 50                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dIndex                 | Distance from 0 to the start of the D region                                                                | 55                    | 54                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| n2Index                | Distance from 0 to the start of the N2 region                                                               | 57                    | 60                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jIndex                 | Distance from 0 to the start of the J region                                                                | 60                    | 62                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| estimatedNumberGenomes | Estimated number of T/B cell genomes present in the sample                                                  | 13                    | 2                     |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| sequenceStatus         | Whether the nucleotide sequence generates a functional amino acid sequence                                  | In                    | In                    |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| cloneResolved          | The gene segments used to construct the CDR3 (VDJ, VJ, DJ)                                                  | VDJ                   | VDJ                   |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vOrphon                | Not implemented at this time                                                                                |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dOrphon                | Not implemented at this time                                                                                |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jOrphon                | Not implemented at this time                                                                                |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| vFunction              | Not implemented at this time                                                                                |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| dFunction              | Not implemented at this time                                                                                |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+
| jFunction              | Not implemented at this time                                                                                |                       |                       |
+------------------------+-------------------------------------------------------------------------------------------------------------+-----------------------+-----------------------+


VDJDB Data Files
################





TCRparse
########

Describe how to parse raw fastq data using legacy tcrdist tools.


Clones File
###########

Mappers
#######

The typical workflow for getting data into tcrdist2 is to load a data file
with Pandas and then apply a mapper to select and rename columns to match
those required by tcrdist2. see :ref:`Mappers`


We have begun writing wrappers for various data sources. Mappers are
just python dictionaries and can be easily emulated if you data headers
are unique.

.. code-block:: python

  import pandas as pd
  from tcrdist import mappers
  mapping = mappers.tcrdist_clone_df_to_tcrdist2_mapping
  tcrdist1_df = mappers.generic_pandas_mapper(df = tcrdist_clone_df_PA,
                                              mapping = mapping)


Required Fields
###############

The follow are required fields for tcrdist2 modules:
