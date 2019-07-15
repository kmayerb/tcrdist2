"""
Contains external functions readPairedSequences(), computeProbs(), and identifyClones()
and internal function for processing TCR sequences invoked by them.

External Functions
------------------

- 'computeProbs' -- Compute probabilities for each TCR based on observed  likelihood from an unpaired dataset
    and based on rearrangement likelihood.
- 'identifyClones' -- Finds unique clones among a paired sequence dataset
- 'readPairedSequences' -- !!!NEEDS FULL DOC!!! - Read a TSV of paired-chain TCR data.

Internal Functions
------------------
- 'samplerProb' -- Compute the probability for the CDR3 protein/nucleotide sequence.
- 'processNT' -- Process one nucleotide TCR sequence (any chain)
- 'getRepresentativeGenes' -- WRAPPER NEEDS DOC -  return util.get_rep(v_gene, organism)
- 'getTCRID' --  Generates a string that is unique for the TCR useful for indexing pairwise distance calculations
    and pd.DataFrames like psDf and clonesDf
"""

# external package imports
import pandas as pd
import numpy as np
from functools import partial
# from .blast import parse_unpaired_dna_sequence_blastn, get_qualstring
# from .objects import TCRChain, TCRClone

    # INTERNAL NOTE
    # for readability tcrdist2 embraces
    # from tcrdist import module syntax
    # inline calls to these functions are module.function(), if we want to change that later we can.
    # examples:
    # tcr_sampler.tcr_sampler()
    # translation.get_translation()
    # compute_probs.rearrangementProb()
    # find_clones.findClones()

from tcrdist import objects
from tcrdist import blast
from tcrdist import util
from tcrdist import tcr_sampler   # added to include tcr_sampler in the namespace
from tcrdist import translation   # added to include get_translation in the namespace
from tcrdist import compute_probs # added to include rearrangementProb in the namespace
from tcrdist import find_clones   # added to include findClones in the namespace

# from .compute_probs import *   # This causes an fatal error, so turned off for intial testing
# from .find_clones import findClones

__all__ = ['processNT',
           'readPairedSequences',
           'filterOutRow',
           'samplerProb',
           'getMaskedSeqs',
           'rearrangementProb',
           'computeProbs',
           'identifyClones',
           'getRepresentativeGenes',
           'getTCRID']


def readPairedSequences(organism, paired_seqs_file, use_parasail = True, try_parasail = True):
    """Read a TSV of paired-chain TCR data.
    Can also accommodate unpaired data or cells with more than two chains
    (e.g. gamma1 and gamma2)

    Parameters
    ----------
    organism : string "human" or "mouse"
    paired_seqs_file (MUST BE TAB DELIMITED .tsv flat file)

    Returns
    -------
    psDf : pd.DataFrame with columns described in Notes

    Notes
    -----

         'id',
         'epitope',
         'subject',
         'a_good_hits',
         'a_status',
         'cdr3a',
         'cdr3a_nucseq',
         'cdr3a_plus',
         'cdr3a_quals',
         'ja_alignlen',
         'ja_bitscore_gap',
         'ja_blast_hits',
         'ja_countreps',
         'ja_evalue',
         'ja_gene',
         'ja_genes',
         'ja_mismatches',
         'ja_mm',
         'ja_rep',
         'ja_reps',
         'va_alignlen',
         'va_bitscore_gap',
         'va_blast_hits',
         'va_countreps',
         'va_evalue',
         'va_gene',
         'va_genes',
         'va_mismatches',
         'va_mm',
         'va_rep',
         'va_reps',
         'b_good_hits',
         'b_status',
         'cdr3b',
         'cdr3b_nucseq',
         'cdr3b_plus',
         'cdr3b_quals',
         'jb_alignlen',
         'jb_bitscore_gap',
         'jb_blast_hits',
         'jb_countreps',
         'jb_evalue',
         'jb_gene',
         'jb_genes',
         'jb_mismatches',
         'jb_mm',
         'jb_rep',
         'jb_reps',
         'vb_alignlen',
         'vb_bitscore_gap',
         'vb_blast_hits',
         'vb_countreps',
         'vb_evalue',
         'vb_gene',
         'vb_genes',
         'vb_mismatches',
         'vb_mm',
         'vb_rep',
         'vb_reps',
         'organism',
         'TCRID

    """
    if use_parasail:
        print("RESULTS BASED ON PARASAIL\n")
    else:
        print("RESULTS BASED ON BLAST")

    raw = pd.read_csv(paired_seqs_file, delimiter='\t')
    """E.g. ['a', 'b']"""
    chains = [s.split('_')[0] for s in raw.columns if s.find('nucseq') > 0]

    out = []
    for c in chains:
        out.append(
            raw.apply(lambda row: processNT(organism=organism,
                                            chain=c.upper(),
                                            nuc=row['%s_nucseq' % c],
                                            quals=row['%s_quals' % c],
                                            use_parasail=use_parasail,
                                            try_parasail=try_parasail).to_series(),
                      axis=1))

    otherCols = [c for c in raw.columns if c.find('nucseq') == -1 and c.find('quals') == -1]
    out = [raw[otherCols]] + out
    psDf = pd.concat(out, axis=1)
    psDf.loc[:, 'organism'] = organism
    psDf.loc[:, 'TCRID'] = psDf.apply(partial(getTCRID, chains=''.join(chains).upper()), axis=1)
    return psDf


def computeProbs(psDf, add_masked_seqs=True, filterOut=False, max_cdr3_length=30, allow_stop_codons=False,
                 allow_X=False):
    """Compute probabilities for each TCR based on observed likelihood from an unpaired dataset
    and based on rearrangement likelihood.

    The resulting pd.DataFrame shares an index with psDf and can be joined:

    psDf = psDf.join(probDf)


    Parameters
    ----------
    psDf : pd.DataFrame
        Paired-sequences dataset generated by td.processing.readPairedSequences()

    Returns
    -------
    probsDf : pd.DataFrame
        Input dataset with columns: a_protseq_prob, cdr3a_protseq_prob, va_rep_prob,
                                    ja_rep_prob, a_nucseq_prob, b_protseq_prob,
                                    cdr3b_protseq_prob, vb_rep_prob, jb_rep_prob, b_nucseq_prob

        optionally includes columns: cdr3a_protseq_masked, a_indels, cdr3a_new_nucseq,
                                     cdr3b_protseq_masked, b_indels, cdr3b_new_nucseq"""

    out = []
    for rowi, row in psDf.iterrows():
        """If iterrows is slow there are potentially ways to speed this up using psDf.apply()"""
        vals = {}
        vals['ind'] = rowi

        if filterOut:
            fo = filterOutRow(row,
                              max_cdr3_length=max_cdr3_length,
                              allow_stop_codons=allow_stop_codons,
                              allow_X=allow_X)
            if fo:
                """vals will be missing keys, which will be assigned Nan in outDf"""
                continue

        aprob_nucseq, aprob_protseq = samplerProb(row, 'a')
        va_rep_prob, ja_rep_prob = compute_probs.rearrangementProb(row, 'a')

        vals['a_protseq_prob'] = aprob_protseq * va_rep_prob * ja_rep_prob
        vals['cdr3a_protseq_prob'] = aprob_protseq
        vals['va_rep_prob'] = va_rep_prob
        vals['ja_rep_prob'] = ja_rep_prob
        vals['a_nucseq_prob'] = aprob_nucseq * va_rep_prob * ja_rep_prob

        bprob_nucseq, bprob_protseq = samplerProb(row, 'b')
        vb_rep_prob, jb_rep_prob = compute_probs.rearrangementProb(row, 'b')

        vals['b_protseq_prob'] = bprob_protseq * vb_rep_prob * jb_rep_prob
        vals['cdr3b_protseq_prob'] = bprob_protseq
        vals['vb_rep_prob'] = vb_rep_prob
        vals['jb_rep_prob'] = jb_rep_prob
        vals['b_nucseq_prob'] = bprob_nucseq * vb_rep_prob * jb_rep_prob

        if add_masked_seqs:
            cdr3a_protseq_masked, ita, cdr3a_new_nucseq = compute_probs.getMaskedSeqs(row, 'a')
            cdr3b_protseq_masked, itb, cdr3b_new_nucseq = compute_probs.getMaskedSeqs(row, 'b')

            vals['cdr3a_protseq_masked'] = cdr3a_protseq_masked
            vals['a_indels'] = ita
            vals['cdr3a_new_nucseq'] = cdr3a_new_nucseq
            vals['cdr3b_protseq_masked'] = cdr3b_protseq_masked
            vals['b_indels'] = itb
            vals['cdr3b_new_nucseq'] = cdr3b_new_nucseq
        out.append(vals)

    outDf = pd.DataFrame(out).set_index('ind')
    assert outDf.shape[0] == psDf.shape[0]
    return outDf


## -- ## copied directly from github (d9394fa9b7)
def samplerProb(r, chain):
    """Compute the probability for the CDR3 protein/nucleotide sequence.

    Parameter
    ---------
    r : pd.Series or TCRChain or TCRClone
        Row from the psDf or any other representation of a TCR chain
    chain : str
        Value 'a' or 'b' or 'A' or 'B'

    Returns
    -------
    prob_nucseq : float
        Probability based on the nucleotide sequence
    prob_protseq : float
        Probability based on the protein sequence"""

    c = chain.lower()

    if c == 'a':
        samplerFunc = tcr_sampler.alpha_cdr3_protseq_probability
    elif c == 'b':
        samplerFunc = tcr_sampler.beta_cdr3_protseq_probability

    prob_nucseq, new_cdr3_nucseq = samplerFunc(r.id,
                                               r.organism,
                                               r['v%s_gene' % c],
                                               r['j%s_gene' % c],
                                               cdr3_protseq='',
                                               cdr3_nucseq=r['cdr3%s_nucseq' % c],
                                               return_final_cdr3_nucseq=True)

    if new_cdr3_nucseq != r['cdr3%s_nucseq' % c]:  ## note note note
        logger.warning('new_cdr3%s_nucseq: %s %s', c, len(new_cdr3_nucseq), new_cdr3_nucseq)
        logger.warning('old_cdr3%s_nucseq: %s %s', c, len(r['cdr3%s_nucseq' % c]), r['cdr3%s_nucseq' % c])
        new_cdr3_protseq = translation.get_translation(new_cdr3_nucseq, '+1')[0]
    else:
        new_cdr3_protseq = r['cdr3%s' % c]
        assert new_cdr3_protseq == translation.get_translation(r['cdr3%s_nucseq' % c], '+1')[0]

    prob_protseq = samplerFunc(r.id,
                               r.organism,
                               r['v%s_gene' % c],
                               r['j%s_gene' % c],
                               cdr3_protseq=new_cdr3_protseq)
    return prob_nucseq, prob_protseq
## -- ## copied directly from github (d9394fa9b7)


def processNT(organism, chain, nuc, quals, use_parasail = True, try_parasail = True):
    """Process one nucleotide TCR sequence (any chain).

    Parameters
    ----------
    organism : str
    chain : str
        Values: alpha, beta, gamma, delta
    nuc : str
    quals : str
        Dot-separated list of integer scores

    Returns
    -------
    TCRChain object
    """

    ch = chain.lower()

    quals = np.array(quals.split('.')).astype(int)
    res = blast.parse_unpaired_dna_sequence_blastn(organism = organism,
                                                   ab = chain,
                                                   blast_seq = nuc,
                                                   info='',
                                                   try_parasail = try_parasail,
                                                   use_parasail = use_parasail,
                                                   nocleanup=True,
                                                   hide_nucseq=False,
                                                   extended_cdr3=True,
                                                   return_all_good_hits=True,
                                                   max_bit_score_delta_for_good_hits=50)
    genes, evalues, status, all_good_hits_with_scores = res
    labels = ['v%s_gene', 'v%s_rep', 'v%s_mm', 'j%s_gene', 'j%s_rep', 'j%s_mm', 'cdr3%s_plus']
    tmp = {g:v for g, v in zip([lab % ch for lab in labels], genes)}
    tmp.update({'%s_evalue' % k.lower():evalues[k][0] for k in list(evalues.keys())})
    tmp.update({'%s_bitscore_gap' % k.lower():evalues[k][1] for k in list(evalues.keys())})

    tmp['%s_status' % ch] = 'OK' if not status else status
    tmp['%s_good_hits' % ch] = all_good_hits_with_scores

    tmp['cdr3%s' % ch], tmp['cdr3%s_nucseq' % ch] = tmp['cdr3%s_plus' % ch].split('-')
    tmp['cdr3%s_quals' % ch] = blast.get_qualstring( tmp['cdr3%s_plus' % ch], nuc, quals )
    tmp['v%s_mismatches' % ch] = tmp['v%s_mm' % ch][0]
    tmp['j%s_mismatches' % ch] = tmp['j%s_mm' % ch][0]
    tmp['v%s_alignlen' % ch] = np.sum(tmp['v%s_mm' % ch])
    tmp['j%s_alignlen' % ch] = np.sum(tmp['j%s_mm' % ch])

    hits = tmp['%s_good_hits' % ch]
    if  hits and len(hits) == 2 and hits[0] and hits[1]:
        tmp['v%s_blast_hits' % ch] = ';'.join( '{}:{}'.format(x[0], x[1]) for x in hits[0] )
        tmp['j%s_blast_hits' % ch] = ';'.join( '{}:{}'.format(x[0], x[1]) for x in hits[1] )
        va_genes = util.get_top_genes( tmp['v%s_blast_hits' % ch] ) ## a set
        ja_genes = util.get_top_genes( tmp['j%s_blast_hits' % ch] )
        tmp['v%s_genes' % ch] = ';'.join( sorted( va_genes ) )
        tmp['j%s_genes' % ch] = ';'.join( sorted( ja_genes ) )
        tmp['v%s_reps' % ch]  = ';'.join( sorted( util.get_top_reps( tmp['v%s_blast_hits' % ch], organism ) ) )
        tmp['j%s_reps' % ch]  = ';'.join( sorted( util.get_top_reps( tmp['j%s_blast_hits' % ch], organism ) ) )
        tmp['v%s_countreps' % ch] = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x, organism) for x in va_genes ))))
        tmp['j%s_countreps' % ch] = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x, organism) for x in ja_genes ))))

    chain = objects.TCRChain(**tmp)
    return chain

def identifyClones(psDf, min_quality_for_singletons=20, average_clone_scores=[], none_score_for_averaging=9.6):
    """

    Finds unique clones among a paired sequence dataset.

    Parameters
    ----------
    psDf : pd.DataFrame
        Paired-sequence dataset generated by readPairedSequences

    min_quality_for_singletons : int

    average_clone_scores : ???? looks like a list

    none_score_for_averaging : float


    Returns
    -------
    clonesDf : pd.DataFrame
        Dataset of (effectively) unique clones in psDf, therefore length
        is <= that of psDf and they do not share an index.

    """

    clonesDf = find_clones.findClones(psDf,
                          min_quality_for_singletons=min_quality_for_singletons,
                          average_clone_scores=average_clone_scores,
                          none_score_for_averaging=none_score_for_averaging)
    chains = [s.split('_')[0] for s in psDf.columns if 'nucseq' in s]
    clonesDf.loc[:, 'CLONEID'] = clonesDf.apply(partial(getTCRID, chains=''.join(chains).upper()), axis=1)
    return clonesDf

def getRepresentativeGenes(organism, v_gene):
    return util.get_rep(v_gene, organism)

def getTCRID(tcr, chains):
    """Generates a string that is unique for the TCR
    useful for indexing pairwise distance calculations
    and pd.DataFrames like psDf and clonesDf

    Parameters
    ----------
    tcr : pd.Series or td.objects.TCRClone
    chains : str
        E.g. 'A' or 'AB' or 'GD'

    Returns
    -------
    s : str
        A string that is unique to the specified TCR"""

    chainSpecificCols = ['v%s_reps', 'cdr3%s', 'j%s_reps']

    cols = [c % chain.lower() for chain in chains for c in chainSpecificCols] + ['organism']

    s = '|'.join([tcr[c] for c in cols if c in tcr])
    return s
