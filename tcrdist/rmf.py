import os
import sys
import pandas as pd

from . import util
from . import parse_tsv
from .all_genes import all_genes
from . import tcr_sampler
from . import paths

def _get_amino_acid_consensus_character( counts ):
    topcount,topaa = max( ( (y,x) for x,y in counts.items() ) )
    frac = float( topcount)/ sum(counts.values())
    if frac >= 0.75:
        return topaa
    elif frac >= 0.4:
        return topaa.lower()
    else:
        return  '.'


def _get_amino_acid_consensus( seqs ):
    numseqs = len(seqs)
    L = len(seqs[0])
    consensus = ''
    for i in range(L):
        counts = {}
        for seq in seqs:
            assert len(seq) == L
            counts[seq[i]] = counts.get(seq[i],0)+1
        consensus += get_amino_acid_consensus_character( counts )
    assert len(consensus) == L
    return consensus

def parse_dataframe(df,k = 'epitope'):
    """
    provides each record of a dataframe as a dict

    replaces parse_tsv

    Paramaters
    ----------
    df : DataFrame

    k : str
        column to make an outer dictionary

    Returns
    -------
    df2 : dictionary
        {"epitope": [{record1},{record2}]}

    """
    assert k in df.columns
    assert isinstance(k, str)
    assert isinstance(df, pd.DataFrame)

    epitopes_in_df = list(df[k].unique())
    infos = df.to_dict("record")
    df2 = {y : [x for x in infos if x['epitope'] == y] for y in epitopes_in_df}

    return(df2)

def _generate_tcrs_dict_from_clones_dataframe(df,
                                              epitopes,
                                              organism,
                                              junction_bars=True,
                                              return_as_tuple = True):
    """
    2019-10-16

    Parameters
    ----------
    df : DataFrame
        clone_dataframe
    epitopes : list
        list of epitopes
    junction_bars : bool

    return_as_tuple : bool
        if true, returns tuple for direct unpacking:
            all_tcrs, all_rep2label_rep, all_rep2label_rep_color = \
            generate_read_motif_dict_from_clones_file()
        if false, the three parts are returned as a dictionary
            r = generate_read_motif_dict_from_clones_file()
            x = r['all_tcrs']
            y = r['all_rep2label_rep']
            z = r['all_rep2label_rep_color']
    Returns
    -------
    dict or tuple:
        Containing Three Parts

        Part 1 - 'all_tcrs'
            {"PA" : [elem]}
            elem contains:
            0:'TRAV7-3*01',
            1:'TRAJ33*01',
            2:'TRBV13-1*01',
            3:'TRBJ2-3*01',
            4:'CAVSLDSNYQLIW',
            5:'CASSDFDWGGDAETLYF',
            6:'TRAV7-3*01',
            7:'TRAJ33*01',
            8:'TRBV13-1*01',
            9:'TRBJ2-3*01',
            10:['V',
                'V',
                ...
                'J',
                'J'],
            11: ['V',
                 'V',...
                 'J'])
            12:  {   'clone_id': 'mouse_tcr0072.clone', ....
                     'va_label_rep_color': 'yellow',
                     'ja_label_rep': 'TRAJ33',
                     'ja_label_rep_color': 'blue',
                     'vb_label_rep': 'TRBV13-1',
                     'vb_label_rep_color': 'black',
                     'jb_label_rep': 'TRBJ2-3',
                     'jb_label_rep_color': 'gold'}}

        Part 2 - 'all_rep2label_rep'

            {'PA': {'TRAV7-3*01': 'TRAV7-3',
                  'TRAJ33*01': 'TRAJ33', ...}

        Part 3 - 'all_rep2label_rep_color'
        'all_rep2label_rep_color' : dict
            {'PA': {  'TRAV7-3*01': 'yellow',
                      'TRAJ33*01': 'blue', ...}

    Example
    -------

    Notes
    -----
    This replaces parse_tsv method in
    _generate_read_motif_dicts_from_clones_file
    """
    all_tcr_infos = parse_dataframe(df, k = "epitope")
    all_tcrs = {}
    all_rep2label_rep = {}
    all_rep2label_rep_color = {}
    for epitope in epitopes:
        infos = all_tcr_infos[epitope]
        util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( infos, organism )

        all_tcrs[epitope] = []
        all_rep2label_rep[epitope] = {}
        all_rep2label_rep_color[epitope] = {}

        for l in infos:

            epitope = l['epitope']
            va_gene = l['va_gene']
            ja_gene = l['ja_gene']
            vb_gene = l['vb_gene']
            jb_gene = l['jb_gene']
            va = l['va_rep']
            ja = l['ja_rep']
            vb = l['vb_rep']
            jb = l['jb_rep']
            cdr3a = l['cdr3a']
            cdr3b = l['cdr3b']
            cdr3a_nucseq = l['cdr3a_nucseq']
            cdr3b_nucseq = l['cdr3b_nucseq']

            ## note-- we are using mm1 reps here, same as in find_cdr3_motifs.py
            ##
            va_rep = all_genes[organism][va].mm1_rep
            ja_rep = all_genes[organism][ja].rep
            vb_rep = all_genes[organism][vb].mm1_rep
            jb_rep = all_genes[organism][jb].rep

            a_junction_results = tcr_sampler.analyze_junction( organism, va_gene, ja_gene, cdr3a, cdr3a_nucseq,
                                                               return_cdr3_nucseq_src=True )
            b_junction_results = tcr_sampler.analyze_junction( organism, vb_gene, jb_gene, cdr3b, cdr3b_nucseq,
                                                               return_cdr3_nucseq_src=True )

            cdr3a_new_nucseq, cdr3a_protseq_masked, cdr3a_protseq_new_nucleotide_countstring,\
                a_trims,a_inserts,cdr3a_nucseq_src = a_junction_results
            cdr3b_new_nucseq, cdr3b_protseq_masked, cdr3b_protseq_new_nucleotide_countstring,\
                b_trims,b_inserts,cdr3b_nucseq_src = b_junction_results

            assert len(cdr3a_nucseq_src) == 3*len(cdr3a)
            assert len(cdr3b_nucseq_src) == 3*len(cdr3b)

            assert type(cdr3b_nucseq_src) == type([])## note

            if junction_bars: ## try to distinguish between N before D and N after D
                for i in range(len(cdr3b_nucseq_src)):
                    if cdr3b_nucseq_src[i] == 'N':
                        if cdr3b_nucseq_src[:i].count('D')==0:
                            cdr3b_nucseq_src[i] = 'N1'
                        else:
                            cdr3b_nucseq_src[i] = 'N2'

            ## let's reconstruct where everything comes from in the cdr3a and cdr3b sequences
            all_tcrs[ epitope ].append( ( va, ja, vb, jb, cdr3a, cdr3b, va_rep, ja_rep, vb_rep, jb_rep,
                                          cdr3a_nucseq_src, cdr3b_nucseq_src, l ) )

            all_rep2label_rep[ epitope ][ va_rep ] = l['va_label_rep']
            all_rep2label_rep[ epitope ][ ja_rep ] = l['ja_label_rep']
            all_rep2label_rep[ epitope ][ vb_rep ] = l['vb_label_rep']
            all_rep2label_rep[ epitope ][ jb_rep ] = l['jb_label_rep']
            all_rep2label_rep_color[ epitope ][ va_rep ] = l['va_label_rep_color']
            all_rep2label_rep_color[ epitope ][ ja_rep ] = l['ja_label_rep_color']
            all_rep2label_rep_color[ epitope ][ vb_rep ] = l['vb_label_rep_color']
            all_rep2label_rep_color[ epitope ][ jb_rep ] = l['jb_label_rep_color']

    if return_as_tuple:
        # return as a tuple for direct unpacking
        return all_tcrs, all_rep2label_rep, all_rep2label_rep_color
    else:
        # return as a dictionary packed
        return({'all_tcrs': all_tcrs,
                'all_rep2label_rep' : all_rep2label_rep,
                'all_rep2label_rep_color' :all_rep2label_rep_color})


def _generate_read_motif_dicts_from_clones_file(clones_file,
                                                epitopes,
                                                organism,
                                                junction_bars=True,
                                                return_as_tuple = False):
    """
    2019-10-16

    Parameters
    ----------
    clones_file : str
        full clones_file for all epitopes in epitopes param
    epitopes : list
        list of epitopes
    junction_bars : bool

    return_as_tuple : bool
        if true, returns tuple for direct unpacking:
            all_tcrs, all_rep2label_rep, all_rep2label_rep_color = \
            generate_read_motif_dict_from_clones_file()
        if false, the three parts are returned as a dictionary
            r = generate_read_motif_dict_from_clones_file()
            x = r['all_tcrs']
            y = r['all_rep2label_rep']
            z = r['all_rep2label_rep_color']
    Returns
    -------
    dict or tuple:
        Containing Three Parts

        Part 1 - 'all_tcrs'
            {"PA" : [elem]}
            elem contains:
            0:'TRAV7-3*01',
            1:'TRAJ33*01',
            2:'TRBV13-1*01',
            3:'TRBJ2-3*01',
            4:'CAVSLDSNYQLIW',
            5:'CASSDFDWGGDAETLYF',
            6:'TRAV7-3*01',
            7:'TRAJ33*01',
            8:'TRBV13-1*01',
            9:'TRBJ2-3*01',
            10:['V',
                'V',
                ...
                'J',
                'J'],
            11: ['V',
                 'V',...
                 'J'])
            12:  {   'clone_id': 'mouse_tcr0072.clone', ....
                     'va_label_rep_color': 'yellow',
                     'ja_label_rep': 'TRAJ33',
                     'ja_label_rep_color': 'blue',
                     'vb_label_rep': 'TRBV13-1',
                     'vb_label_rep_color': 'black',
                     'jb_label_rep': 'TRBJ2-3',
                     'jb_label_rep_color': 'gold'}}

        Part 2 - 'all_rep2label_rep'

            {'PA': {'TRAV7-3*01': 'TRAV7-3',
                  'TRAJ33*01': 'TRAJ33', ...}

        Part 3 - 'all_rep2label_rep_color'
        'all_rep2label_rep_color' : dict
            {'PA': {  'TRAV7-3*01': 'yellow',
                      'TRAJ33*01': 'blue', ...}

        Example
        -------

    """
    # all_tcr_infos - THIS IS EQUIVALENT TO THE CLONES FILE AS A DICTIONARY KEYED ON EPITOPE : COLUMN
    all_tcr_infos = parse_tsv.parse_tsv_file( clones_file, ['epitope'], [], True )
    all_tcrs = {}
    all_rep2label_rep = {}
    all_rep2label_rep_color = {}
    for epitope in epitopes:
        infos = all_tcr_infos[epitope]
        util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( infos, organism )

        all_tcrs[epitope] = []
        all_rep2label_rep[epitope] = {}
        all_rep2label_rep_color[epitope] = {}

        for l in infos:

            epitope = l['epitope']
            va_gene = l['va_gene']
            ja_gene = l['ja_gene']
            vb_gene = l['vb_gene']
            jb_gene = l['jb_gene']
            va = l['va_rep']
            ja = l['ja_rep']
            vb = l['vb_rep']
            jb = l['jb_rep']
            cdr3a = l['cdr3a']
            cdr3b = l['cdr3b']
            cdr3a_nucseq = l['cdr3a_nucseq']
            cdr3b_nucseq = l['cdr3b_nucseq']

            ## note-- we are using mm1 reps here, same as in find_cdr3_motifs.py
            ##
            va_rep = all_genes[organism][va].mm1_rep
            ja_rep = all_genes[organism][ja].rep
            vb_rep = all_genes[organism][vb].mm1_rep
            jb_rep = all_genes[organism][jb].rep

            a_junction_results = tcr_sampler.analyze_junction( organism, va_gene, ja_gene, cdr3a, cdr3a_nucseq,
                                                               return_cdr3_nucseq_src=True )
            b_junction_results = tcr_sampler.analyze_junction( organism, vb_gene, jb_gene, cdr3b, cdr3b_nucseq,
                                                               return_cdr3_nucseq_src=True )

            cdr3a_new_nucseq, cdr3a_protseq_masked, cdr3a_protseq_new_nucleotide_countstring,\
                a_trims,a_inserts,cdr3a_nucseq_src = a_junction_results
            cdr3b_new_nucseq, cdr3b_protseq_masked, cdr3b_protseq_new_nucleotide_countstring,\
                b_trims,b_inserts,cdr3b_nucseq_src = b_junction_results

            assert len(cdr3a_nucseq_src) == 3*len(cdr3a)
            assert len(cdr3b_nucseq_src) == 3*len(cdr3b)

            assert type(cdr3b_nucseq_src) == type([])## note

            if junction_bars: ## try to distinguish between N before D and N after D
                for i in range(len(cdr3b_nucseq_src)):
                    if cdr3b_nucseq_src[i] == 'N':
                        if cdr3b_nucseq_src[:i].count('D')==0:
                            cdr3b_nucseq_src[i] = 'N1'
                        else:
                            cdr3b_nucseq_src[i] = 'N2'

            ## let's reconstruct where everything comes from in the cdr3a and cdr3b sequences
            all_tcrs[ epitope ].append( ( va, ja, vb, jb, cdr3a, cdr3b, va_rep, ja_rep, vb_rep, jb_rep,
                                          cdr3a_nucseq_src, cdr3b_nucseq_src, l ) )

            all_rep2label_rep[ epitope ][ va_rep ] = l['va_label_rep']
            all_rep2label_rep[ epitope ][ ja_rep ] = l['ja_label_rep']
            all_rep2label_rep[ epitope ][ vb_rep ] = l['vb_label_rep']
            all_rep2label_rep[ epitope ][ jb_rep ] = l['jb_label_rep']
            all_rep2label_rep_color[ epitope ][ va_rep ] = l['va_label_rep_color']
            all_rep2label_rep_color[ epitope ][ ja_rep ] = l['ja_label_rep_color']
            all_rep2label_rep_color[ epitope ][ vb_rep ] = l['vb_label_rep_color']
            all_rep2label_rep_color[ epitope ][ jb_rep ] = l['jb_label_rep_color']

    if return_as_tuple:
        # return as a tuple for direct unpacking
        return all_tcrs, all_rep2label_rep, all_rep2label_rep_color
    else:
        # return as a dictionary packed
        return({'all_tcrs': all_tcrs,
                'all_rep2label_rep' : all_rep2label_rep,
                'all_rep2label_rep_color' :all_rep2label_rep_color})


def _generate_read_motif_ng_tcrs_dict(chains,
                                     organism = "mouse",
                                     ng_log_file_A = None,
                                     ng_log_file_B = None,
                                     ng_log_file_G = None,
                                     ng_log_file_D = None,
                                     ng_path       = None,
                                     max_ng_lines  = 5000000):
    """
    Parameters
    ----------
    chains : list
        ["A", "B"] or ["G", "D"] or ["A"]
    organism : str
        "mouse" or "human"
    ng_log_file_A : str
        optional name of ng_reference file ALPHA chain TCRs
    ng_log_file_B : str
        optional name of ng_reference file BETA chain TCRs
    ng_log_file_G : str
        optional name of ng_reference file GAMMA chain TCRs
    ng_log_file_D : str
        optional name of ng_reference file DELTA chain TCRs
    ng_path : str
        optional specify location of ng_reference

    Returns
    -------
    ng_tcrs : dict


    """
    # If no ng_log_file_X are specified, then use defaults
    if ng_log_file_A is None:
        ng_log_file_A = 'new_nextgen_chains_{}_A.tsv'.format(organism)
    if ng_log_file_B is None:
        ng_log_file_B = 'new_nextgen_chains_{}_B.tsv'.format(organism)
    if ng_log_file_G is None:
        ng_log_file_G = 'new_nextgen_chains_{}_G.tsv'.format(organism)
    if ng_log_file_D is None:
        ng_log_file_D = 'new_nextgen_chains_{}_D.tsv'.format(organism)
    if ng_path is None:
        ng_path = paths.path_to_current_db_files()

    ng_tcrs = {k:[] for k in chains}

    ## index these by the v_rep and the j_rep
    for chain in chains:
        if chain == "A":
            ng_logfile = ng_log_file_A
        if chain == "B":
            ng_logfile = ng_log_file_B
        if chain == "G":
            ng_logfile = ng_log_file_G
        if chain == "D":
            ng_logfile = ng_log_file_D

        ng_logfile = os.path.join(ng_path, ng_logfile)

        if not os.path.isfile(ng_logfile):
            raise OSError('WARNING::missing nextgen TCR chains file: {}'.format(ng_logfile))

        counter=0
        num_chains=0
        ab_chains = {}

        for line in open(ng_logfile,'r'):
            counter+=1
            if max_ng_lines and counter>max_ng_lines:break
            l = line[:-1].split('\t')
            if counter==1:
                assert l==['v_reps','j_reps','cdr3','cdr3_nucseq']
                continue
            if not counter % 500000:
                pass
                #basic.Log(`counter`+' '+`num_chains`+' '+ng_logfile)
            v_reps = set( ( util.get_mm1_rep( x, organism ) for x in l[0].split(',') ) ) ## mm1 reps
            j_reps = l[1].split(',')
            cdr3,cdr3_nucseq = l[2:4]

            ## now add to the different places
            for v_rep in v_reps:
                for j_rep in j_reps:
                    if v_rep not in ab_chains: ab_chains[v_rep] = {}
                    if j_rep not in ab_chains[v_rep]: ab_chains[v_rep][j_rep] = []
                    ab_chains[v_rep][j_rep].append( (cdr3, cdr3_nucseq ))

            num_chains += 1
        sys.stdout.write('read {} {}-chains from {}\n'.format(num_chains,chain,ng_logfile))
        ng_tcrs[chain] = ab_chains
    return(ng_tcrs)

def create_wtd_pwm_from_sequences( seqs, alphabet, target_reps, reps ):
    """
    seqs :

    alpahabet :

    target_reps :

    reps :

    """
    assert len(seqs) == len(reps)
    num_target_reps = len(target_reps)
    ## now iteratively adjust the wts on rep-pairs
    for bigrepeat in range(10):
        reppair_wts = {}
        if bigrepeat==0:
            for rp in reps:
                reppair_wts[rp] = 1.0 ## starting guess
        else:
            for rp in reps:
                reppair_wts[rp] = 0.75 + 0.5 * random.random() ## starting guess

        prev_dev = 1e6
        for repeat in range(100):

            ## what's the deviation
            dev = 0.0
            for ii in range(2):
                ii_target_reps = [x[ii] for x in target_reps]
                ii_reps = [x[ii] for x in reps]

                scale_factor = float( len(reps ) )/ len(target_reps)

                counts = {}
                for rp in reps:
                    counts[rp[ii]] = counts.get(rp[ii],0) + reppair_wts[rp]

                for rep,count in counts.items():
                    desired_count = scale_factor * ii_target_reps.count(rep)
                    dev += abs( desired_count - count )
                    fac = float(desired_count)/count
                    adjust = fac**(0.25)

                    #print 'desired_count:',desired_count,'count:',count,'fac:',fac,'adjust:',adjust,rep

                    for rp in reppair_wts:
                        if rp[ii] == rep:
                            reppair_wts[rp] *= adjust
            #print 'repeat:',repeat,'dev:',dev
            if abs(prev_dev-dev)<1e-3 and dev<1e-1:
            #if abs(dev)<1e-1:
                break
            prev_dev = dev

        #print 'final_dev:', bigrepeat,dev
        if dev<1e-1:
            break


    L = len(seqs[0])
    pwm = {}
    for i in range(L):
        pwm[i] = dict(zip(alphabet,[0.0]*len(alphabet)))

    for seq,rp in zip( seqs, reps ):
        assert len(seq) == L
        seqwt = reppair_wts[ rp ]
        #print seq, rp, seqwt
        for i,a in enumerate(seq):
            pwm[i][a] += seqwt

    for i in range(L):
        tot = sum( pwm[i].values() )
        for a in alphabet:
            pwm[i][a] /= tot
    return pwm


def analyze_matches_using_ngseqs( matches, matched_tcrs, ab, tcrs):
    #global tcrs
    import numpy as np
    ng_lenseqs = []
    ng_fwdseqs = []
    ng_revseqs = []

    ng_fwdseq_reps = []
    ng_lenseq_reps = []
    matched_reps = []

    seen = set() ## no repeats of ngseqs
    seen_samelen = set() ## no repeats of ngseqs

    for (mseq,nseq,positions,rpositions),ii in zip( matches, matched_tcrs ):
        tcr = tcrs[ii]
        if ab == 'A':
            my_cdr3,vrep,jrep = tcr[4:5]+tcr[6: 8]
        else:
            my_cdr3,vrep,jrep = tcr[5:6]+tcr[8:10]
        matched_reps.append( ( vrep, jrep ) )

        mylen = len(my_cdr3)
        if vrep in ng_tcrs[ab] and jrep in ng_tcrs[ab][vrep]:
            ngl = [ x for x in ng_tcrs[ab][vrep][jrep] if x not in seen ]
            if not ngl:
                print('empty ngl!')

            for ngseq in random.sample( ngl, min(num_nextgen_samples,len(ngl)) ):
                seen.add(ngseq)
                (cdr3,cdr3_nucseq) = ngseq
                L = len(cdr3)
                fseq = ''
                rseq = ''
                for pos in positions:
                    if pos>=L:
                        fseq += '-'
                    else:
                        fseq += cdr3[pos]
                for pos in rpositions:
                    if pos>=L:
                        rseq += '-'
                    else:
                        rseq += cdr3[L-1-pos]
                ng_fwdseqs.append(fseq)
                ng_revseqs.append(rseq)
                ng_fwdseq_reps.append( ( vrep, jrep ) )

            ## cdr3s with the same length
            ngl_samelen = [ x for x in ng_tcrs[ab][vrep][jrep] if len(x[0]) == mylen and x not in seen_samelen ]
            if not ngl_samelen:
                print('empty ngl_samelen!')
            for ngseq in random.sample( ngl_samelen, min(num_nextgen_samples,len(ngl_samelen))):
                seen_samelen.add( ngseq )
                cdr3 = ngseq[0]
                ng_lenseqs.append( ''.join( [ cdr3[x] for x in positions ] ) )
                ng_lenseq_reps.append( ( vrep, jrep ) )

    pwm = logo_tools.create_protein_pwm_from_sequences( [x[0] for x in matches ])

    npwm_alphabet = junction_bars_order[ab] if junction_bars else ['V','N','D','J']
    npwm = logo_tools.create_pwm_from_sequences( [x[1] for x in matches ], npwm_alphabet )

    #nbr_pwm = logo_tools.create_protein_pwm_from_sequences( [x[0] for x in nbr_matches ])
    #nbr_npwm = logo_tools.create_pwm_from_sequences( [x[1] for x in nbr_matches ], ['V','D','J','N'] )

    if ng_lenseqs:
        ng_lenpwm = create_wtd_pwm_from_sequences( ng_lenseqs, amino_acids+['-'], matched_reps, ng_lenseq_reps )
    else:
        ng_lenpwm = 0

    ng_fwdpwm = create_wtd_pwm_from_sequences( ng_fwdseqs, amino_acids+['-'], matched_reps, ng_fwdseq_reps )
    ng_revpwm = create_wtd_pwm_from_sequences( ng_revseqs, amino_acids+['-'], matched_reps, ng_fwdseq_reps )

    N = len(pwm)
    fwdpwm = {}
    revpwm = {}
    for i in range(N):
        fwdpwm[i] = {}
        revpwm[i] = {}
        incrememnt = 1.0/len(matches)
        for pos in [x[2][i] for x in matches]:
            fwdpwm[i]['pos'] = fwdpwm[i].get('pos',0)+incrememnt
        for pos in [x[3][i] for x in matches]:
            revpwm[i]['pos'] = revpwm[i].get('pos',0)+incrememnt

    ## look at relative entropies between nbrpwm and the fwd and rev pwms
    ## not nbr anymore since this is a subroutine
    ##
    scale_by_relent = {}
    for i in range(N):
        relents=[]
        for control_pwm in [ ng_fwdpwm[i], ng_revpwm[i] ]:
            relent = 0.0
            for a,pa in pwm[i].items():
                if pa>= min_prob_for_relent_for_scaling:
                    qa = max(min_prob_for_relent_for_scaling, control_pwm.get(a,min_prob_for_relent_for_scaling))
                    paqa = np.log2(pa/qa)
                    relent += pa * paqa
            relents.append( relent )
        scale_by_relent[i] = max(0.,min(1., min(relents)/max_column_relent_for_scaling) )
        print('RE {:2d} {:5.2f} {:5.2f} {:5.2f} {} {} {}'.format( i, min(relents), relents[0], relents[1], ab, epitope, ''.join(showmotif) ))

    return pwm, npwm, ng_lenpwm, ng_fwdpwm, ng_revpwm, fwdpwm, revpwm, scale_by_relent, \
        ng_fwdseq_reps, ng_lenseq_reps, len( ng_lenseqs ), len( ng_fwdseqs )

## make a v-gene logo
def get_counts_list_condensing_alleles( counts_string, rep2label_rep, rep2label_rep_color ):
    #global rep2label_rep
    #global rep2label_rep_color
    counts ={}
    for tag,count in [x.split(':') for x in counts_string.split(',') ]:
        rc = ( rep2label_rep[ tag ][4:], rep2label_rep_color[ tag ] )
        counts[rc] = counts.get(rc,0)+float(count)
    return [ (y,x[0],x[1]) for x,y in counts.items() ]

def get_counts_lists_from_tcr_indices( indices ):
    vcounts = {}
    jcounts = {}
    for ii in indices:
        tcr = tcrs[ii]
        if ab == 'A':
            vrep,jrep = tcr[6: 8]
        else:
            vrep,jrep = tcr[8:10]
        vcounts[vrep] = vcounts.get(vrep,0)+1
        jcounts[jrep] = jcounts.get(jrep,0)+1
    vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.items()] )
    jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.items()] )
    return get_counts_list_condensing_alleles(vstring, rep2label_rep, rep2label_rep_color),\
            get_counts_list_condensing_alleles(jstring, rep2label_rep, rep2label_rep_color)

def get_counts_lists_from_rep_lists( reps ):
    vcounts = {}
    jcounts = {}
    for vrep,jrep in reps:
        vcounts[vrep] = vcounts.get(vrep,0)+1
        jcounts[jrep] = jcounts.get(jrep,0)+1
    vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.items()] )
    jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.items()] )
    return get_counts_list_condensing_alleles(vstring, rep2label_rep, rep2label_rep_color),\
            get_counts_list_condensing_alleles(jstring, rep2label_rep, rep2label_rep_color)


## KMB: THE IMPORTANT OUTPUT OF THIS BLOCK SEEMS TO BE < all_neighbors >
def generate_all_neighbors(clones_file,
                           chains = ["A", "B"],
                           nbr_distance = 100,
                           epitope = None,
                           dist_file_A = None,
                           dist_file_B = None,
                           dist_file_G = None,
                           dist_file_D = None):
    """

    Returns
    -------

    all_neighbors : dict
        dict of lists


    """
    if "A" in chains and dist_file_A is None:
        dist_file_A = '{}_A_{}.dist'.format(clones_file[:-4], epitope)
    if "B" in chains and  dist_file_B is None:
        dist_file_B = '{}_B_{}.dist'.format(clones_file[:-4], epitope)
    if "G" in chains and  dist_file_G is None:
        dist_file_G = '{}_G_{}.dist'.format(clones_file[:-4], epitope)
    if "D" in chains and  dist_file_D is None:
        dist_file_D = '{}_D_{}.dist'.format(clones_file[:-4], epitope)

    all_distances = {}
    all_neighbors = {}
    for ab in chains:
        if ab is 'A':
            distfile = dist_file_A
        if ab is 'B':
            distfile = dist_file_B
        if ab is 'G':
            distfile = dist_file_G
        if ab is 'D':
            distfile = dist_file_D

        assert os.path.isfile(distfile)
        N=0
        all_nbrs = []
        all_dists = []
        for line in open( distfile,'r'):
            l = line.split()
            clone_id = l[0]
            index = len(all_dists)
            #assert tcrs[ index ][-1]['clone_id'] == clone_id
            dists = [ float(x) for x in l[1:] ]
            if not N:
                N = len(dists)
            else:
                assert N == len(dists)

            nbrs = []
            for ii,d in enumerate(dists):
                if d <= nbr_distance:
                    nbrs.append( ii )
            all_dists.append( dists )
            all_nbrs.append( nbrs )

        #assert len(all_nbrs) == len(tcrs)

        all_distances[ab] = all_dists
        all_neighbors[ab] = all_nbrs
    return(all_neighbors)

    #
def generate_all_nbr_from_dataframe(dist_df, nbr_distance = 100.0):
    """
    Given a distance matrix, return a list of index positions of all neighbors
    within a certain distance.



    Parameters
    ----------
    dist_df : DataFrame

    nbr_distanace : float

    Returns
    -------
    all_nbrs : list
        list of lists

    Examples
    --------
    >>> M = pd.DataFrame([[0,2,5,1],
    ...                   [2,0,3,1],
    ...                   [5,3,0,5],
    ...                   [1,1,5,0]])
    >>> rmf.generate_all_nbr_from_dataframe(M,2)
    [0, 1, 3], [0, 1, 3], [2], [0, 1, 3]]


    Notes
    ----
    Methods for iterating a DataFrame like a flat file
    replace:
        for line in open( distfile,'r'):
    with:
    option 1
    x = distA.to_dict("record")
    for r in x:
        print(list(r.values()))

    option 2
    for index, row in distA.iterrows():
        print(index, row.to_list())
    """
    N=0
    all_nbrs = []
    all_dists = []
    for ind, row in dist_df.iterrows():
        l = row.to_list()
        clone_id = ind
        index = len(all_dists)
        #assert tcrs[ index ][-1]['clone_id'] == clone_id
        dists = [ float(x) for x in l ]
        if not N:
            N = len(dists)
        else:
            assert N == len(dists)

        nbrs = []
        for ii,d in enumerate(dists):
            if d <= nbr_distance:
                nbrs.append( ii )
        all_dists.append( dists )
        all_nbrs.append( nbrs )

    return(all_nbrs)


## make a v-gene logo
def get_counts_list_condensing_alleles( counts_string, rep2label_rep, rep2label_rep_color ):
    #global rep2label_rep
    #global rep2label_rep_color
    counts ={}
    for tag,count in [x.split(':') for x in counts_string.split(',') ]:
        rc = ( rep2label_rep[ tag ][4:], rep2label_rep_color[ tag ] )
        counts[rc] = counts.get(rc,0)+float(count)
    return [ (y,x[0],x[1]) for x,y in counts.items() ]

def get_counts_lists_from_tcr_indices( indices, tcrs, chain ):
    vcounts = {}
    jcounts = {}
    for ii in indices:
        tcr = tcrs[ii]
        if chain in ['A', "G"]:
            vrep,jrep = tcr[6: 8]
        else:
            vrep,jrep = tcr[8:10]
        vcounts[vrep] = vcounts.get(vrep,0)+1
        jcounts[jrep] = jcounts.get(jrep,0)+1
    vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.items()] )
    jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.items()] )
    return get_counts_list_condensing_alleles(vstring, rep2label_rep, rep2label_rep_color),\
            get_counts_list_condensing_alleles(jstring, rep2label_rep, rep2label_rep_color)

def get_counts_lists_from_rep_lists( reps ):
    vcounts = {}
    jcounts = {}
    for vrep,jrep in reps:
        vcounts[vrep] = vcounts.get(vrep,0)+1
        jcounts[jrep] = jcounts.get(jrep,0)+1
    vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.items()] )
    jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.items()] )
    return get_counts_list_condensing_alleles(vstring, rep2label_rep, rep2label_rep_color),\
            get_counts_list_condensing_alleles(jstring, rep2label_rep, rep2label_rep_color)
