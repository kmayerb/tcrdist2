"""
This page contains documentation for functions used for preprocessing
Adaptive Biotechnology's Immunoseq Files. Currently only single
chain (Beta) data is supported by these functions.

Examples
--------

1. Preprocess Bulk Beta-Chain Dataset

.. code-block:: python

    from tcrdist.all_genes import all_genes
    from tcrdist import preprocess_adaptive
    from tcrdist.repertoire import TCRrep

    adaptive_df = preprocess_adaptive.parse_immunoseq(organism = 'human', filename = "example1.tsv" )
    nonnull_cdr3 = adaptive_df.cdr3_b_aa != ""
    adaptive_df = adaptive_df[nonnull_cdr3].copy()

    tr = TCRrep(adaptive_df, organism = "human", chains = ['beta'])
    tr.infer_cdrs_from_v_gene(chain = 'beta', imgt_aligned=True)
    tr.index_col = tr.index_cols = ['id', 'subject',
                   'v_b_gene', 'j_b_gene',
                   'cdr3_b_aa','cdr3_b_nucseq',
                   'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']


2. Get a record of gene names in the adaptive file relative to the IMGT
assigned name.

.. code-block:: python

    from tcrdist import preprocess_adaptive
    v_b_genes, j_b_genes, v_b_genes_imgt, j_b_genes_imgt = preprocess_adaptive.\
    _get_adaptive_gene_names("example1.tsv", organism = "human")

"""


from collections import namedtuple, Counter
import numpy as np
import pandas as pd
from .all_genes import all_genes

def parse_immunoseq(filename, organism, subject = None):
    """
    Creates a cell_df DataFrame compatible with tcrdist2 using
    an Adaptive Biotechnologies immunoSEQ file as the source.

    Parameters
    ----------
    filename : str
        filename of .tab deliminited Adaptive Biotechnologies immunoSEQ file
        with standard headers.
    orgnanism : str
        'human' or 'mouse'
    subject : str
        The subject associated with this file

    Returns
    -------
    immunoseq_cell_df : DataFrame

    Notes
    -----
    !!!! WARNING: when the Adpative Gene call does not specify an allele
    the *01 allele is assigned.

    See page 24 section 3.4.13 SAMPLE EXPORT of the immunoSEQ instrument
     `Analyzer Manual <https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf>`_

    This currently only parses the beta-chain

    Some confusing nomenclature in tcrdist

    * v_b_gene is the beta-chain V gene at the highest possible level of resolution, including allele if available
    * vb_rep is the allele used to represent that gene family (e.g. *02 and *01 both are repressented as by *01)
    * vb_countreps is the gene level, ignoring allele, representative used for counting/
    """
    if subject is None:
        subject = str(filename.replace(".tsv", ""))

    storage_list = list()
    with open(filename , "r") as fh:
        # First Line, We check the header against the namedtuple
        _validate_header_of_adaptive_immunoseq_file(header = fh.readline(),
                                                    adaptive = adaptive)
        # Subsequent Lines
        for i, line in enumerate(fh):
            n = len(adaptive._fields)
            # only take the first n elements
            elements = line.strip().split("\t")[0:n]
            # initialize adaptive named tuple
            adapt = adaptive(*elements)

            # get v_b_gene, vb_countreps names
            v_b_gene, vb_countreps = _get_v_gene_attrs(adapt)
            j_b_gene, jb_countreps = _get_j_gene_attrs(adapt)

            # _
            # Adaptive's Gene Names Need to Match all_genes (IMGT Name)
            #    TCRBV05-07*01 --- >  TRBV5-7*01
            # The steps are:
            # 1.  TCRB -> TRBC (e.g., TCRBV05-07*01 ->  TRBV05-7*01 )
            # 2.  e.g., TRBV05-07*01 ->  TRBV05-7*01
            # 3.  TRBV05-7*01 -> TRBV5-7*01
            # 4.  If no allele is called, then we assign it '*01'
            # 5.  There are some '-1' names not recognized
            # (e.g., TRBV19-1*01 to TRBV19*01)
            # Replace Adaptive's TCRB naming convention in favor of IMGT name
            v_b_gene = v_b_gene.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BV0","BV")
            j_b_gene = j_b_gene.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BJ0","BJ")

            vb_countreps = vb_countreps.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BV0","BV")
            jb_countreps = jb_countreps.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BV0","BV")

            # Assign allele *01 to genes not resolved at the allele level
            if v_b_gene.find('*0') == -1:
                v_b_gene = v_b_gene + "*01"
            if j_b_gene.find('*0') == -1:
                j_b_gene = j_b_gene + "*01"

            # TRBV19-1*01 to TRBV19*01
            if v_b_gene not in all_genes[organism].keys():
                v_b_gene = v_b_gene.replace("-1","")

            if j_b_gene not in all_genes[organism].keys():
                j_b_gene = j_b_gene.replace("-1","")

            # get the primary amino acid sequence sequence
            cdr3_b_aa = adapt.aminoAcid
            cdr3_b_nucseq = _extract_cdr3_nucseq_from_adaptive(adapt)


            storage_list.append( [i,
                                  subject,
                                  v_b_gene,
                                  vb_countreps,
                                  j_b_gene,
                                  jb_countreps,
                                  cdr3_b_aa,
                                  cdr3_b_nucseq ])

        immunoseq_cell_df = pd.DataFrame(storage_list, columns = ['id',
                                                                  'subject',
                                                                  'v_b_gene',
                                                                  'vb_countreps',
                                                                  'j_b_gene',
                                                                  'jb_countreps',
                                                                  'cdr3_b_aa',
                                                                  'cdr3_b_nucseq'])
        return immunoseq_cell_df

"""
These are the fields we expect to find in and Adaptive Biotechnologies
immunoSEQ file. (See end of this file for a complete table)

See page 24 section 3.4.13 SAMPLE EXPORT of the immunoSEQ instrument
 `Analyzer Manual <https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf>`_
"""
adaptive = namedtuple("adaptive_sample",
                      [
                        "nucleotide",
                        "aminoAcid",
                        "count",
                        "frequencyCount",
                        "cdr3Length",
                        "vMaxResolved",
                        "vFamilyName",
                        "vGeneName",
                        "vGeneAllele",
                        "vFamilyTies",
                        "vGeneNameTies",
                        "vGeneAlleleTies",
                        "dMaxResolved",
                        "dFamilyName",
                        "dGeneName",
                        "dGeneAllele",
                        "dFamilyTies",
                        "dGeneNameTies",
                        "dGeneAlleleTies",
                        "jMaxResolved",
                        "jFamilyName",
                        "jGeneName",
                        "jGeneAllele",
                        "jFamilyTies",
                        "jGeneNameTies",
                        "jGeneAlleleTies",
                        "vDeletion",
                        "n1Insertion",
                        "d5Deletion",
                        "d3Deletion",
                        "n2Insertion",
                        "jDeletion",
                        "vIndex",
                        "n1Index",
                        "dIndex",
                        "n2Index",
                        "jIndex",
                        "estimatedNumberGenomes",
                        "sequenceStatus",
                        "cloneResolved",
                        "vOrphon",
                        "dOrphon",
                        "jOrphon",
                        "vFunction",
                        "dFunction",
                        "jFunction"])


def _validate_header_of_adaptive_immunoseq_file(header, adaptive):
    """
    Checks that the element of the header of input file
    match the expected field names in a numedtuple.

    Parameters
    ----------
    header : str
        tab-delimited string representing the first header line
    adaptive : collections.namedtuple
        named tuple of fields expected in an adaptive file

    Returns
    -------
    boolean
        if header fields are validit

    Raises
    ------
    ValueError
        if any of the header fields do not match fields in expected adaptive file

    Notes
    -----
    See the 3.4.13 SAMPLE EXPORT section of immunoSEQ
    `Analyzer Manual <https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf>`_


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

    """
    # Check validity of file header against expected fields
    n_fields = len(adaptive._fields)

    header = header.strip().split("\t")
    header = [x.split(" ")[0] for x in header]
    # only take the portion of the header names before the first space

    compare_header_to_adaptive = \
        [(i,j) for i,j in zip(header[0:n_fields], adaptive._fields) if i!=j]

    if len(compare_header_to_adaptive) > 0:
        print("ERROR: THESE FIELDS DO NOT MATCH EXPECTATION")
        print(compare_header_to_adaptive)

    if not np.all([i==j for i,j in zip(header[0:n_fields], adaptive._fields) ] ):
        raise ValueError("YOUR FILE DOES NOT MATCH THE EXPECTED FIRST " +
                         "{} COLUMNS IN AN ADAPTIVE IMMUNOSEQ FILE".format(n_fields))

    return 1


def _extract_cdr3_nucseq_from_adaptive(adapt):
    """
    Extracts the nucleotide sequence coding the CDR3, defined as starting at
    canonical V-gene Cystine Residue

    Parameters
    ----------
    adapt : collections.namedtuple
        instance of named tuple of fields expected in an adaptive file

    Returns
    -------
    cdr3_b_nucseq : str

    Notes
    -----

    ADAPTIVE's output does not provide explict indices for the cdr3 nucleotide
    sequence. Thus, we use the length of the amino acid sequence
    and the start index of cys_residue (see page 18 of their
    manual for 'v_index' defining the start of Cys residue, 'jIndex' is
    start of J gene.

    """
    length_cdr3_nucseq = len(adapt.aminoAcid)*3
    cys_residue    = int(adapt.vIndex)
    start_J_gene   = int(adapt.jIndex)

    length_partial_nuc = len(adapt.nucleotide[cys_residue:start_J_gene])
    remaining_length_nuc = length_cdr3_nucseq - length_partial_nuc
    end_index = start_J_gene + remaining_length_nuc

    cdr3_nucseq = adapt.nucleotide[cys_residue:end_index]
    return cdr3_nucseq


def _get_v_gene_attrs(adapt):
    """
    Returns v-gene attributes: Allele and Gene Name.

    Parameters
    ----------
    adapt : collections.namedtuple
        instance of named tuple of fields expected in an adaptive file

    Returns
    -------
    v_b_gene : str
        gene at highest possible level of resolution, including allele if available
    vb_countreps : str
        the count representative used is the gene level, ignoring allele.
    """
    if adapt.vGeneName is not "":
        v_b_gene     = adapt.vMaxResolved
        vb_countreps = adapt.vGeneName
    # if no unique vGeneName was assigned, check ties, and take the first
    # Adaptive output sorts them in ascending order (i.e. TCRBV06-01,TCRBV06-05,TCRBV06-06)
    elif adapt.vGeneNameTies is not "":
        v_b_gene = adapt.vGeneNameTies.split(",")[0]
        vb_countreps   = adapt.vGeneNameTies.split(",")[0]
    else:
        v_b_gene = None
        vb_countreps = None
    return (v_b_gene, vb_countreps)


def _get_j_gene_attrs(adapt):
    """
    Returns j-gene attributes: Allele and Gene Name.

    Parameters
    ----------
    adapt : collections.namedtuple
        instance of named tuple of fields expected in an adaptive file

    Returns
    -------
    vj_b_gene : str
        gene at highest possible level of resolution, including allele if available
    jb_countreps : str
        is the gene level, ignoring allele, representative used

    Notes
    -----
    Assign j_b_gene similar to above v_b gene
    """
    if adapt.jGeneName is not "":
        j_b_gene = adapt.jMaxResolved
        jb_countreps  = adapt.jGeneName
    elif adapt.jGeneNameTies is not "":
        j_b_gene = adapt.jGeneNameTies.split(",")[0]
        jb_countreps = adapt.jGeneNameTies.split(",")[0]
    else:
        j_b_gene     = None
        jb_countreps = None
    return (j_b_gene, jb_countreps)

def _get_adaptive_gene_names(filename, organism):
    """"
    Gets raw adaptive gene_names from an immunoseq file.

    Parameters
    ----------
    filenames : Adaptive SAMPLE EXPORT of the immunoSEQ instrument

    organism : str
        'human' or 'mouse'

    Return
    ------
    Counter(v_b_genes_raw) : collections.Counter
    Counter(j_b_genes_raw) : collections.Counter
    v_b_genes_imgt : dict
    j_b_genes_imgt : dict
    """
    v_b_genes_raw = list()
    j_b_genes_raw = list()
    v_b_genes_imgt = dict()
    j_b_genes_imgt = dict()

    with open(filename , "r") as fh:
        # First Line, We check the header against the namedtuple
        _validate_header_of_adaptive_immunoseq_file(header = fh.readline(),
                                                    adaptive = adaptive)
        # Subsequent Lines
        for i, line in enumerate(fh):
            n = len(adaptive._fields)
            # only take the first n elements
            elements = line.strip().split("\t")[0:n]
            # initialize adaptive named tuple
            adapt = adaptive(*elements)
            v_b_gene, vb_countreps = _get_v_gene_attrs(adapt)
            j_b_gene, jb_countreps = _get_j_gene_attrs(adapt)

            # Save the Raw Genes
            v_b_gene_raw = str(v_b_gene)
            j_b_gene_raw = str(j_b_gene)

            v_b_genes_raw.append(v_b_gene_raw)
            j_b_genes_raw.append(j_b_gene_raw)

            v_b_gene = v_b_gene.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BV0","BV")
            j_b_gene = j_b_gene.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BJ0","BJ")

            vb_countreps = vb_countreps.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BV0","BV")
            jb_countreps = jb_countreps.replace("TCRB", "TRB").\
                replace("-0", "-").replace("BV0","BV")

            # Assign allele *01 to genes not resolved at the allele level
            if v_b_gene.find('*0') == -1:
                v_b_gene = v_b_gene + "*01"
            if j_b_gene.find('*0') == -1:
                j_b_gene = j_b_gene + "*01"

            # TRBV19-1*01 to TRBV19*01
            if v_b_gene not in all_genes[organism].keys():
                v_b_gene = v_b_gene.replace("-1","")

            if j_b_gene not in all_genes[organism].keys():
                j_b_gene = j_b_gene.replace("-1","")

            v_b_genes_imgt[v_b_gene_raw] = v_b_gene
            j_b_genes_imgt[j_b_gene_raw] = j_b_gene


    return Counter(v_b_genes_raw), Counter(j_b_genes_raw), v_b_genes_imgt, j_b_genes_imgt



"""
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

"""
