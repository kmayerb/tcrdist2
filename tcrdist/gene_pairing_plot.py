import numpy as np
import pandas as pd
import itertools
from scipy import stats

import statsmodels.api as sm
#from sklearn.metrics import adjusted_mutual_info_score

#from basic import *
import html_colors
#import svg_basic

__all__ = ['plotPairings']

greek_alpha = '&#x3b1;'
greek_beta  = '&#x3b2;'

segtype2greek_label = { 'VA':'V'+greek_alpha, 'JA':'J'+greek_alpha,
                        'VB':'V'+greek_beta , 'JB':'J'+greek_beta }

def _create_file(cmds, width, height, filename, background_color=None, use_xlink=False):
    out = open(filename, 'w')
    extra = '' if not use_xlink else 'xmlns:xlink="http://www.w3.org/1999/xlink"'
    out.write('<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg" version="1.1" {} >\n'\
              .format(int(width),int(height),extra))
    if background_color:
        out.write(_rectangle( (0,0), (width,height), background_color, 'white', 0 ) )
    out.write('\n'.join(cmds)+'\n')
    out.write('</svg>\n')
    out.close()

def _rectangle( upper_left, lower_right, fill, stroke, stroke_width=1, dashed=False ):
    stroke_dasharray_style= ';stroke-dasharray:5,5' if dashed else ''
    return '<rect x="{}" y="{}" height="{}" width="{}" style="fill:{};stroke:{};stroke-width:{}{}" />\n'\
        .format( upper_left[0], upper_left[1], lower_right[1]-upper_left[1], lower_right[0]-upper_left[0],
                 fill, stroke, stroke_width, stroke_dasharray_style )

def _make_text( text, lower_left, fontsize,
               extra_tag = None, font_family = "monospace", color = "black", font_weight = "normal" ):
    assert font_weight in ['normal','bold']
    cmd = '<text x="{:.3f}" y="{:.3f}" font-size="{}" font-weight="{}" font-family="{}" fill="{}" xml:space="preserve">{}</text>\n'\
        .format( lower_left[0], lower_left[1], fontsize, font_weight, font_family, color, text )
    return cmd

def _linear_gradient_cmd(gradient_id_counter, x1, y1, x2, y2, offsets, colors, spread_method="pad" ):
    defs_id_prefix = 'TCRGRAD'
    gradient_id_counter += 1
    id = '{}lingrad{:d}'.format(defs_id_prefix, gradient_id_counter)

    stoplines = ''
    assert len(offsets) == len(colors)

    for offset,color in zip( offsets,colors):
        stoplines += """<stop offset="{:.1f}%"   stop-color="{}" stop-opacity="1"/>""".format( offset, color )

    cmd = """<defs>
                    <linearGradient id="{}"
                    x1="{:.1f}%" y1="{:.1f}%"
                    x2="{:.1f}%" y2="{:.1f}%"
                    spreadMethod="{}">
                                    {}
                    </linearGradient>
            </defs>""".format( id, x1, y1, x2, y2, spread_method, stoplines )

    return gradient_id_counter, id, cmd

def _computeAssociations(df, cols, countCol='Count'):
    res = []
    col1, col2 = cols
    for val1, val2 in itertools.product(df[col1].unique(), df[col2].unique()):
        OR, pvalue, tab = _testAssociation(df, (col1, val1), (col2, val2), countCol=countCol)
        tot = np.sum(tab)
        res.append({'OR':OR,
                    'pvalue':pvalue,
                    'Col0':col1,
                    'Col1':col2,
                    'Val0':val1,
                    'Val1':val2,
                    '00':tab[0, 0] / tot,
                    '01':tab[0, 1] / tot,
                    '10':tab[1, 0] / tot,
                    '11':tab[1, 1] / tot})
    resDf = pd.DataFrame(res)
    resDf.loc[:, 'qvalue'] = sm.stats.multipletests(resDf['pvalue'].values, method='fdr_bh')[1]
    resDf = resDf.sort_values(by='pvalue', ascending=True)
    return resDf

def _testAssociation(df, node1, node2, countCol='Count'):
    """Test if the occurence of nodeA paired with nodeB is more/less common than expected.

    Parameters
    ----------
    nodeX : tuple (column, value)
        Specify the node by its column name and the value.

    Returns
    -------
    OR : float
        Odds-ratio associated with the 2x2 contingency table
    pvalue : float
        P-value associated with the Fisher's exact test that H0: OR = 1"""
    col1, val1 = node1
    col2, val2 = node2
    
    tmp = df[[col1, col2, countCol]].dropna()
    #print(node1, node2, countCol)
    #print(tmp)

    tab = np.zeros((2, 2))
    tab[0, 0] = (((tmp[col1]!=val1) & (tmp[col2]!=val2)) * tmp[countCol]).sum()
    tab[0, 1] = (((tmp[col1]!=val1) & (tmp[col2]==val2)) * tmp[countCol]).sum()
    tab[1, 0] = (((tmp[col1]==val1) & (tmp[col2]!=val2)) * tmp[countCol]).sum()
    tab[1, 1] = (((tmp[col1]==val1) & (tmp[col2]==val2)) * tmp[countCol]).sum()

    """Add 1 to cells with zero"""
    tab[tab==0] = 1
    #OR, pvalue = fisherTest(tab)
    OR, pvalue = stats.fisher_exact(tab)
    return OR, pvalue, tab

def plotPairings(filename, df, cols, countCol=None, use_color_gradients=True, other_frequency_threshold=0.01):
    """Diagram depicts the gene-segment pairing structure of the dataset. The four
    genes (cols) are arrayed left to right. Below each gene-type label (eg "VA")
    is a color-stack showing all the TCR clones and how they break down into the different genes for that gene-type. Each clone
    is devoted a constant vertical height in pixels indicated in the text at the top (N pixels in "Nx y-pixel scale"). The curved
    segments joining neighboring gene-stacks show how the two gene distributions pair up, with the thickness of the segments
    corresponding to the number of clones having those two segments (scaled by the indicated y-pixel scale). Significant gene-gene
    pairings (positive or negative correlations with a P-value less than 1e-6) are labeled at the beginning and ending of the
    corresponding segments. Gene-gene pairings which are not observed and for which this under-representation is significant
    are indicated by dashed segments with P-value labels. Enrichments (depletions) of gene segments relative to
    background are shown for all labeled genes by up (down) arrows where the number of arrowheads reflects the base-2
    logarithm of the fold change, rounded down (one arrowhead means 2 <= fold change < 4,
    two arrowheads means 4 <= fold change < 8, and so on).

    The left-right ordering of the segment types is chosen so that VA and JA are on the left, VB and JB are on the right,
    and the alpha-beta pairing with the largest adjusted mutual information is in the middle.

    Parameters
    ----------
    filename : str
        Name of SVG file for storing the output figure.
    df : pd.DataFrame
        Contains gene segment data, one row per clone.
    cols : list
        List of columns for displaying frequency and pairings, in order from left to right.
    countCol : str
        Optionally provide a count or frequency column for weighting the rows/clones

    Returns
    -------
    filename : str
        File name if output is successfully created."""

    params = dict(min_ami_for_colorscale=0.114, # from shuffling experiments
                    max_ami_for_colorscale=0.5,
                    min_entropy_for_colorscale=0.0,
                    max_entropy_for_colorscale=5.0,
                    min_jsd_for_colorscale=0.02259, ## from background vs background comparisons
                    max_jsd_for_colorscale=0.0,
                    min_gene_frequency_for_labels=0.05,
                    reverse_gradients=False,
                    pval_threshold_for_plotting_gene_correlations=1e-2,
                    pval_threshold_for_svg_correlations=1e-6)

    gradient_id_counter = 0

    """Pixel units"""
    left_margin = 50
    right_margin = 50
    top_margin = 50
    bottom_margin = 50
    yspacer = 50

    flat_band = 50
    final_flat_band = flat_band if use_color_gradients else 2.5*flat_band
    middle_band = 400
    slope_weight = 100
    pairing_svg_y_offset = top_margin

    ff = 'Droid Sans Mono'

    # ypixel_scale = max(1, int( 0.5 + 600.0/df.shape[0] ) )
    ypixel_scale = 100
    pairing_svg_width = left_margin + right_margin + 3*(flat_band+middle_band) + final_flat_band

    epitope_fontsize = 60
    midpoint = left_margin + 2*flat_band + 1.5*middle_band
    pairing_svg_cmds = []
    '''
    pairing_svg_cmds.append( _make_text( '{}'.format(epitope, len(tcrs) ),
                                                  [midpoint-0.5*0.6*epitope_fontsize*len(epitope),
                                                   pairing_svg_y_offset+epitope_fontsize-20], epitope_fontsize,
                                                  font_family=ff ) )
    '''

    correlation_fontsize = 14.
    correlation_fontheight = correlation_fontsize * 0.75

    if countCol is None:
        df = df.assign(Count=1)
        countCol = 'Count'
    
    """Compute marginal frequencies for each value within each column"""
    tot = df[countCol].sum()
    freqs = {}
    for c in cols:
        """Assign "Other" to all values below a threshold"""
        tmp = df.groupby(c)[countCol].agg(lambda v: np.sum(v) / tot).sort_values(ascending=False)
        ind = tmp.index[tmp < other_frequency_threshold]
        df.loc[df[c].isin(ind), c] = 'Other'
        freqs[c] = df.groupby(c)[countCol].agg(lambda v: np.sum(v) / tot).sort_values(ascending=False)
        print(c)
        print(freqs[c])

    """For each consecutive pair of columns compute 2 x 2 table of frequencies, OR and Fisher's Exact p-value"""
    res = []
    for i in range(len(cols) - 1):
        res.append(_computeAssociations(df, cols=(cols[i], cols[i+1]), countCol=countCol))
    resDf = pd.concat(res, axis=0, ignore_index=True)
    resDf = resDf.sort_values(by='pvalue').set_index(['Val0', 'Val1'])
    print(resDf)

    for ii in range(len(cols) - 1):
        correlation_paths = []
        r0 = cols[ii]
        r1 = cols[ii+1]

        x0 = left_margin + ii*( flat_band + middle_band )

        text = segtype2greek_label[ r0 ]
        fontsize = 40.
        xtext = x0 + 0.5*flat_band - 0.5*0.6*fontsize*2
        ytext = pairing_svg_y_offset+yspacer-6
        ## hacking
        ytext -= 6
        if ii==0:
            xtext += 8
        pairing_svg_cmds.append( _make_text( text, [ xtext, ytext ], fontsize, font_family=ff ) )
        if ii == (len(cols) - 2): ## add the final column label
            text = segtype2greek_label[ r1 ]
            xtext = x0+1.5*flat_band-0.5*0.6*fontsize*2+middle_band
            xtext -= 8
            pairing_svg_cmds.append( _make_text( text, [ xtext, ytext ], fontsize, font_family=ff ) )
        
        '''
        pairing_svg_cmds.append( _make_text( '(AMI: {:.2f})'.format(ami),
                                                      [x0+flat_band+middle_band/2.5, pairing_svg_y_offset+yspacer-5],
                                                      12, font_family=ff  ))
        '''

        
        #vl = [ (y,x) for x,y in repcounts[r0].iteritems() ]
        #jl = [ (y,x) for x,y in repcounts[r1].iteritems() ]

        #vl.sort() ; vl.reverse()
        #jl.sort() ; jl.reverse()

        #a0colors = dict(zip( [x[1] for x in vl], html_colors.get_rank_colors_no_lights( len(vl) ) ) )
        #a1colors = dict(zip( [x[1] for x in jl], html_colors.get_rank_colors_no_lights( len(jl) ) ) )
        r0colors = dict(zip(freqs[r0].index, html_colors.get_rank_colors_no_lights(len(freqs[r0]))))
        r1colors = dict(zip(freqs[r1].index, html_colors.get_rank_colors_no_lights(len(freqs[r1]))))
        #print(r0colors)
        #print(r1colors)
        '''
        reps2tcrs = {}
        for t in tcrs:
            vj = ( t[ rep_index[r0] ], t[ rep_index[r1] ] )
            if vj not in reps2tcrs:reps2tcrs[vj] = []
            reps2tcrs[vj].append( t )
        '''

        ## on the left, the V-segments, ordered by counts
        ## on the right, J-segments, ordered by counts

        ## need to assign a vertical range to each v/j segment
        ## start with, one pixel per tcr
        ##

        a1counts = {}

        yleft = yspacer + pairing_svg_y_offset
        """Nested loops over the alleles (a0, a1) in the pair of columns, r0, r1"""
        for a0, a0freq in freqs[r0].iteritems():
            y0_right = yspacer + pairing_svg_y_offset
            a0color = r0colors[a0]
            for a1, a1freq in freqs[r1].iteritems():
                a1color = r1colors[a1]
                vj = (a0, a1)
                a1color = r1colors[a1]

                """Frequency of this allele pairing"""
                freq_tcrs = resDf.loc[(a0, a1), '11']
                num_tcrs_scaled = freq_tcrs * ypixel_scale
                stroke_width = np.ceil(num_tcrs_scaled)

                ## ok make a spline
                yright = y0_right + a1counts.get(a1, 0)

                #line/spline points
                j_flat_band = flat_band if ii<(len(cols) - 2) else final_flat_band
                points = [ (np.floor(x0), yleft + 0.5*num_tcrs_scaled ),
                           (x0+flat_band, yleft + 0.5*num_tcrs_scaled ),
                           (np.ceil(x0 + flat_band + middle_band), yright + 0.5*num_tcrs_scaled ),
                           (np.ceil(x0 + flat_band + middle_band + j_flat_band), yright + 0.5*num_tcrs_scaled ) ]


                path1_cmds = 'M {} {} L {} {} M {} {} C {} {}, {} {}, {} {}'\
                    .format( points[0][0], points[0][1], ## start of v-line
                             points[1][0], points[1][1], ## end point of v-line
                             points[1][0], points[1][1],
                             points[1][0] + slope_weight, points[1][1], ## control for spline start
                             points[2][0] - slope_weight, points[2][1], ## control for spline end
                             points[2][0], points[2][1] )

                if use_color_gradients:
                    path1a_cmds = 'M {} {} L {} {}'\
                        .format( points[0][0], points[0][1],  ## start of v-line
                                 points[1][0], points[1][1] ) ## end point of v-line
                    pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                             .format(path1a_cmds, a0color, stroke_width ) )
                    ## define the gradient
                    path1b_cmds = 'M {} {} C {} {}, {} {}, {} {}'\
                        .format( points[1][0], points[1][1],
                                 points[1][0] + slope_weight, points[1][1], ## control for spline start
                                 points[2][0] - slope_weight, points[2][1], ## control for spline end
                                 points[2][0], points[2][1] )
                    #v_line_rhs_fraction = float(flat_band) / (flat_band + middle_band )
                    offsets = [0, 25.0, 75.0, 100]
                    #offsets = [0, 45.0, 55.0, 100]
                    #offsets = [0, 90.0, 99.0, 100]
                    if params['reverse_gradients']:
                        colors = [a1color, a1color, a0color, a0color]
                    else:
                        colors = [a0color, a0color, a1color, a1color]
                    gradient_id_counter, gradient_id, gradient_cmd = _linear_gradient_cmd(gradient_id_counter, 0, 0, 100, 0, offsets, colors )
                    pairing_svg_cmds.append( gradient_cmd )
                    pairing_svg_cmds.append( '<path d="{}" stroke="url(#{})" stroke-width="{}" fill="none"/>'\
                                             .format(path1b_cmds, gradient_id, stroke_width ) )
                else:
                    pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                             .format(path1_cmds,a0color, stroke_width ) )

                if ii == (len(cols) - 2): ## add the right-most flat band
                    path2_cmds = 'M {} {} L {} {}'\
                        .format( points[2][0], points[2][1], ## start of j-line
                                 points[3][0], points[3][1] ) ## end of j-line

                    pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                             .format(path2_cmds, a1color, stroke_width) )

                '''
                if vj in epitope_correlations_svg[epitope] and not paper_figs:
                    #print 'vj has correlations:',vj,epitope_correlations_svg[epitope][vj]
                    if not num_tcrs:
                        #print 'make dotted line!',vj
                        if not ( no_pairing_text or paper_figs ):
                            if paper_supp:
                                assert use_color_gradients
                                ## define the gradient
                                path1b_cmds = 'M {} {} C {} {}, {} {}, {} {}'\
                                    .format( points[1][0], points[1][1],
                                             points[1][0] +slope_weight, points[1][1], ## control for spline start
                                             points[2][0] -slope_weight, points[2][1], ## control for spline end
                                             points[2][0], points[2][1] )
                                #v_line_rhs_fraction = float(flat_band) / (flat_band + middle_band )
                                offsets = [0, 25.0, 75.0, 100]
                                #offsets = [0, 45.0, 55.0, 100]
                                #offsets = [0, 90.0, 99.0, 100]
                                colors = [a0color, a0color, a1color, a1color]
                                gradient_id, gradient_cmd = _linear_gradient_cmd( 0, 0, 100, 0, offsets, colors )
                                pairing_svg_cmds.append( gradient_cmd )
                                pairing_svg_cmds.append( '<path d="{}" stroke="url(#{})" stroke-width="2" stroke-dasharray="5,5" fill="none"/>'\
                                                         .format(path1b_cmds, gradient_id, stroke_width ) )

                            else:
                                dotted_path_cmds = 'M {} {} L {} {} M {} {} C {} {}, {} {}, {} {}'\
                                    .format( points[0][0], points[0][1], ## start of v-line
                                             points[1][0], points[1][1], ## end point of v-line
                                             points[1][0], points[1][1],
                                             points[1][0] +slope_weight, points[1][1], ## control for spline start
                                             points[2][0] -slope_weight, points[2][1], ## control for spline end
                                             points[2][0], points[2][1] )
                                pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="2" stroke-dasharray="5,5" fill="none"/>'\
                                                         .format(path1_cmds,a0color ) )

                    ## new way, just use regular text elements
                    ## pretend that the spline is actually a straight line between these points
                    swf=0.4
                    yshift = correlation_fontheight*0.5
                    p0 = ( points[1][0]+slope_weight*swf, points[1][1]+yshift )
                    p1 = ( points[2][0]-slope_weight*swf, points[2][1]+yshift )

                    dx = p1[0]-p0[0]
                    dy = p1[1]-p0[1]

                    ## so, what is the rotation we'll need?
                    rotangle = math.atan2(dy,dx) * ( 180.0 / math.pi )

                    step = 0.05
                    lower_left  = [ p0[0] + step*dx, p0[1] + step*dy ]

                    step = 0.95
                    lower_right = [ p0[0] + step*dx, p0[1] + step*dy ]

                    pval,enrich = epitope_correlations_svg[epitope][vj]

                    ## write some curved text
                    if enrich==0:
                        msg = '0x ({:.0E})'.format(pval)
                    elif enrich<0.1:
                        msg = '{:.2f}x ({:.0E})'.format(enrich,pval)
                    else:
                        msg = '{:.1f}x ({:.0E})'.format(enrich,pval)

                    fill1,fill2 = 'black','black'
                    if a0color=='black':
                        fill1 = 'gold'
                    if a1color=='black' and use_color_gradients:
                        fill2 = 'gold'


                    cmd1 = '<text x="{:.3f}" y="{:.3f}" font-size="{}" font-family="{}" fill="{}" transform="rotate({:.3f},{:.3f},{:.3f})" >{}</text>\n'\
                        .format( lower_left[0], lower_left[1], correlation_fontsize, ff, fill1,
                                 rotangle, lower_left[0], lower_left[1], msg )

                    cmd2 = '<text text-anchor="end" x="{:.3f}" y="{:.3f}" font-size="{}" font-family="{}" fill="{}" transform="rotate({:.3f},{:.3f},{:.3f})" >{}</text>\n'\
                        .format( lower_right[0], lower_right[1], correlation_fontsize, ff, fill2,
                                 rotangle, lower_right[0], lower_right[1], msg )

                    correlation_paths.append( ( pval, (cmd1, cmd2) ) )
                    #print 'corr cmd1:',vj,cmd1
                    #print 'corr cmd2:',vj,cmd2
                '''

                yleft += num_tcrs_scaled
                a1counts[a1] = a1counts.get(a1, 0) + num_tcrs_scaled
                y0_right += a1freq * ypixel_scale

        '''
        ## try doing the p-val paths
        correlation_paths.sort()
        correlation_paths.reverse() ## go in decreasing order of p-val so the most significant are on top

        ## now write the text
        for (pval, cmds ) in correlation_paths:
            pairing_svg_cmds.extend( cmds )

        ## let's label the alleles in the left stack (and right stack if ii==2)
        fontsize = 40 if paper_figs else 20.0 if paper_supp else 20
        fontheight = 0.75*fontsize
        fontwidth =  0.6 *fontsize

        min_height_for_labels = fontheight+1

        for jj,(r,ll,repcolors) in enumerate( [ (r0,vl,a0colors),(r1,jl,a1colors)]  ):
            if ii<2 and jj>0:continue


            ## label in white?
            x = x0 + jj*(flat_band+middle_band)
            ystart = yspacer+pairing_svg_y_offset

            for ( count,rep) in ll:
                if count*ypixel_scale < min_height_for_labels: break
                #ystop    = ystart + count*ypixel_scale
                midpoint =  ystart + count*ypixel_scale*0.5
                text = rep[2:]

                lower_left = [ x+2, midpoint+fontheight/2.0 ]

                my_flat_band = final_flat_band if ii==2 and jj==1 else flat_band
                bgcolor = repcolors[rep]
                textcolor = 'black' if ((paper_figs or paper_supp) and bgcolor!= 'black') else 'white'
                textcolor = 'black' if bgcolor!= 'black' else 'white'

                if True or paper_figs or paper_supp: ## center the text, unless on either side...
                    text_width = fontwidth*len(text)
                    lower_left_ha = {'left'  : lower_left,
                                     'right' : [ x+my_flat_band-text_width, midpoint+fontheight/2.0 ],
                                     'center': [ x+0.5*my_flat_band-0.5*text_width, midpoint+fontheight/2.0 ]}
                    if jj==0 and ii==0: ## left-most guy
                        ha = 'left'
                    elif jj==1 and ii==2: ## right-most guy
                        ha = 'right'
                    else:
                        ha = 'center'
                    pairing_svg_cmds.append( _make_text( text, lower_left_ha[ha], fontsize, color=textcolor,
                                                                  font_family=ff))


                elif (True or jj==1) and fontwidth*len(text)>my_flat_band: # right-most set, dont want to over-run
                    myfontsize=int(0.5+(my_flat_band-4)/(len(text)*0.6))
                    pairing_svg_cmds.append( _make_text( text, lower_left, myfontsize, color=textcolor,
                                                                  font_family=ff))
                else:
                    pairing_svg_cmds.append( _make_text( text, lower_left, fontsize, color=textcolor,
                                                                  font_family=ff))



                ## add an enrichment glyph?
                if make_enrichment_glyphs:
                    enrich = float( all_countrep_enrichment[ epitope ][ rep ][0][0] )
                    if enrich>=2. or enrich<=0.5:
                        ## add a glyph
                        if paper_supp or paper_figs:
                            arrow_length = 1.35 * min_height_for_labels
                            arrow_width = 3.5
                        else:
                            arrow_length = 1.35 * min_height_for_labels
                            arrow_width = 1.5
                        #arrow_length = min_height_for_labels
                        #arrow_width = 2.5
                        eg_sep = 14.0
                        if 'A' in r:
                            center = [ lower_left_ha[ha][0] + text_width + eg_sep, midpoint ]
                        else:
                            #print rep
                            assert 'B' in r
                            center = [ lower_left_ha[ha][0] - eg_sep, midpoint ]

                        pairing_svg_cmds += svg_basic.enrichment_glyph_cmds( center, arrow_length, arrow_width,
                                                                             enrich )


                ystart += count*ypixel_scale
        '''

        pairing_svg_y_offset += 2*yspacer #+ ypixel_scale

    '''
    if no_pairing_text:
        tmpcmds = pairing_svg_cmds[:]
        pairing_svg_cmds = []
        for cmd in tmpcmds:
            if '<text' in cmd:
                print 'skip:',cmd
            else:
                pairing_svg_cmds.append( cmd )
    '''

    """Make SVG file"""
    bg_color = None # 'white'
    _create_file(pairing_svg_cmds, pairing_svg_width, pairing_svg_y_offset + bottom_margin,
                           filename, background_color=bg_color)
    return filename

if __name__ == '__main__':
    np.random.seed(110820)
    n = 50
    df = pd.DataFrame({'VA':np.random.choice(['TRAV14', 'TRAV12', 'TRAV3', 'TRAV23', 'TRAV11', 'TRAV6'], n),
                       'JA':np.random.choice(['TRAJ4', 'TRAJ2', 'TRAJ3','TRAJ5', 'TRAJ21', 'TRAJ13'], n),
                       'VB':np.random.choice(['TRBV14', 'TRBV12', 'TRBV3', 'TRBV23', 'TRBV11', 'TRBV6'], n),
                       'JB':np.random.choice(['TRBJ4', 'TRBJ2', 'TRBJ3','TRBJ5', 'TRBJ21', 'TRBJ13'], n)})
    df = df.assign(Count=1)
    df.loc[:10, 'Count'] = 10

    print(df)
    plotPairings('test.svg', df, ['JA', 'VA', 'VB'], countCol='Count')

'''
## load epitope jsd values
epitope_jsds = {}
jsd_tsvfile = clones_file[:-4] + '_JS_divergence.tsv'
if not exists( jsd_tsvfile ):
    print 'Sorry, you need to run analyze_gene_frequencies.py before running make_gene_plots.py'
    exit()

lines = parse_tsv_file( jsd_tsvfile, [], ['epitope'] + [x+'_jsd_normed' for x in segtypes_lowercase] )
for line in lines:
    epitope = line[0]
    vals = map(float,line[1:])
    epitope_jsds[epitope] = {}
    assert len(vals)== len(segtypes)
    for segtype,val in zip( segtypes, vals ):
        epitope_jsds[ epitope ][ segtype ] = val

epitope_entropies = {}
epitope_mis = {}
epitope_correlations = {}
epitope_correlations_svg = {}
epitope_repcounts = {}
epitope_repcounts_by_len = {}

all_tcrs = parse_tsv_file( clones_file, ['epitope'], [], True )


## returns id, cmd

pairing_svg_cmds = []
path_def_counter = 0
make_enrichment_glyphs = ( countrep_enrichments_file != None )
if make_enrichment_glyphs:
    all_countrep_enrichment = parse_tsv_file( countrep_enrichments_file, [ 'epitope','gene' ], ['jsd_prob_enrich'] )


if not epitopes:
    epitopes = all_tcrs.keys()[:]
    epitopes.sort()

for epitope in epitopes:

    ## this fills in *_label_rep fields in all_tcrs dictionary
    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( all_tcrs[epitope], organism )

    epitope_entropies[epitope] = {}
    epitope_mis[epitope] = {}
    epitope_correlations[epitope] = []
    epitope_repcounts[epitope] = {}
    epitope_correlations_svg[epitope] = {}

    """tcrs is list of tuples (col1, col2, denom1, denom2)"""
    tcrs = []
    for fulltcr in all_tcrs[epitope]:
        tcrs.append( ( fulltcr['va_label_rep'], fulltcr['ja_label_rep'],
                       fulltcr['vb_label_rep'], fulltcr['jb_label_rep'],
                       len(fulltcr['cdr3a']), len(fulltcr['cdr3b'] ) ) ) # not subtracting 5 any more


    repcounts = {}
    repcounts2 = {}

    repcounts_by_len = {}

    for i,r in enumerate(segtypes):
        repcounts[r] = {}
        repcounts_by_len[r] = {}
        for s in segtypes[i+1:]:
            repcounts2[(r,s)] = {}

    rep_index = dict(zip(segtypes,range(len(segtypes))))

    for tcr in tcrs:
        assert len(tcr) == 6
        for r in segtypes:
            rep = tcr[ rep_index[r] ]
            repcounts[r][rep] = repcounts[r].get(rep,0)+1

            assert r[1] in 'AB'
            cdr3len = tcr[4] if r[1]=='A' else tcr[5]

            min_cdr3len = min(min_cdr3len,cdr3len)
            max_cdr3len = max(max_cdr3len,cdr3len)

            if cdr3len not in repcounts_by_len[r]:
                repcounts_by_len[r][cdr3len] = {}
            repcounts_by_len[r][cdr3len][rep] = repcounts_by_len[r][cdr3len].get(rep,0)+1

        for rs in repcounts2:
            rep = (tcr[ rep_index[rs[0]]], tcr[ rep_index[rs[1]]] )
            repcounts2[rs][rep] = repcounts2[rs].get(rep,0)+1

    for r in segtypes:
        for s in segtypes:
            rs=(r,s)
            if rs in repcounts2:
                for rep1 in repcounts[r]:
                    for rep2 in repcounts[s]:
                        rep=(rep1,rep2)
                        if rep not in repcounts2[rs]:
                            repcounts2[rs][rep] = 0



    epitope_repcounts[epitope] = dict( repcounts )
    epitope_repcounts_by_len[epitope] = dict( repcounts_by_len )

    N = len(tcrs)

    ## compute entropies, mutual informations
    for r in segtypes:
        entropy=0
        for rep,count in repcounts[r].iteritems():
            prob=float(count)/N
            entropy -= prob * math.log(prob,2)
        print 'ENT {:4s} {} entropy: {:7.3f} entropy_pow2: {:7.3f} N: {:6d}'.format(epitope,r,entropy,2**entropy,N)
        epitope_entropies[epitope][r] = entropy


    from sklearn.metrics import adjusted_mutual_info_score
    from scipy.stats import hypergeom

    all_ab_amis = []
    all_amis = {}

    for rs in repcounts2:
        ab_pairing = ( rs[0][1] != rs[1][1] )
        cluster_pairing = ( rs[0][0] == 'C' or rs[1][0] == 'C' )

        mi=0.0
        entropy=0
        for (rep1,rep2),count in repcounts2[rs].iteritems():
            pxy = float(count)/N
            if pxy>0: entropy -= pxy*math.log(pxy,2)
            count1 = repcounts[rs[0]][rep1]
            count2 = repcounts[rs[1]][rep2]
            px = float(count1)/N
            py = float(count2)/N
            if pxy>0: mi += pxy * math.log( (pxy/ (px*py)), 2 )

            ## lets look at the significance of this overlap
            expected = px * py * N
            pval = 1

            if count > expected:
                ## compute hypergeometric distn prob
                max_possible_overlap = min(count1,count2)
                x = np.arange(0,max_possible_overlap+1)
                cdf = hypergeom.cdf( x, N, count1, count2 ) ## cdf is accumulated prob <= val
                sf  = hypergeom.sf( x, N, count1, count2 )
                pval = sf[count-1] ## now greater than or equal to count
                if pval<1e-3:
                    print 'PVAL: {:4s} {:12.3e} {}-{} {:15s} {:15s} overlap: {:4d} expect: {:7.1f} count1: {:4d} count2: {:4d} '\
                        .format(epitope,pval,rs[0],rs[1],str(rep1),str(rep2),count,expected,count1,count2)
                #exit()

                if pval<pval_threshold_for_svg_correlations:
                    #print 'svg pval!',rep1,rep2,pval
                    epitope_correlations_svg[epitope][(rep1,rep2)] = ( pval, count/expected )
                    epitope_correlations_svg[epitope][(rep2,rep1)] = ( pval, count/expected )

            if count < expected:
                ## compute hypergeometric distn prob
                max_possible_overlap = min(count1,count2)
                x = np.arange(0,max_possible_overlap+1)
                cdf = hypergeom.cdf( x, N, count1, count2 ) ## cdf is accumulated prob <= val
                sf  = hypergeom.sf( x, N, count1, count2 )
                pval = cdf[count] ## less than or equal to count
                if pval<1e-3:
                    print 'PVAL: {:4s} {:12.3e} {}-{} {:15s} {:15s} overlap: {:4d} expect: {:7.1f} count1: {:4d} count2: {:4d} '\
                        .format(epitope,pval,rs[0],rs[1],str(rep1),str(rep2),count,expected,count1,count2)
                #exit()
                if pval<pval_threshold_for_svg_correlations:
                    #print 'svg pval!',rep1,rep2,pval
                    epitope_correlations_svg[epitope][(rep1,rep2)] = ( pval, count/expected )
                    epitope_correlations_svg[epitope][(rep2,rep1)] = ( pval, count/expected )


            if ab_pairing and (not cluster_pairing) and pval<pval_threshold_for_plotting_gene_correlations:
                if count==0:
                    logenrich = math.log(  0.25 / expected, 2 )
                else:
                    logenrich = math.log( count / expected, 2 )

                epitope_correlations[epitope].append ( ( logenrich, -1*math.log( pval,10 ), rs, rep1, rep2 ) )


        ## compute an adjusted mutual information score
        labels0 = []
        labels1 = []

        tcr_labels0 = []
        tcr_labels1 = []

        for tcr in tcrs:
            l0 = tcr[ rep_index[ rs[0] ] ]
            l1 = tcr[ rep_index[ rs[1] ] ]
            if l0 not in labels0: labels0.append( l0 )
            if l1 not in labels1: labels1.append( l1 )
            tcr_labels0.append( labels0.index(l0) )
            tcr_labels1.append( labels1.index(l1) )

        ami = adjusted_mutual_info_score( tcr_labels0, tcr_labels1 )


        if ab_pairing:
            all_ab_amis.append( ( ami, rs ) )

        all_amis[ (rs[0],rs[1]) ] = ami
        all_amis[ (rs[1],rs[0]) ] = ami

        print 'MI {:4s} {}-{} MI: {:7.3f} AMI: {:7.3f} MI_pow2 {:7.3f} entropy: {:7.3f} entropy_pow2: {:7.3f}'\
            .format(epitope,rs[0],rs[1],mi,ami,2**mi,entropy,2**entropy)

        epitope_entropies[epitope][rs] = entropy
        epitope_mis[epitope][rs] = (mi,ami)



    all_ab_amis.sort()
    all_ab_amis.reverse()

    top_pairing = all_ab_amis[0]
    print 'top ab pairing:',top_pairing

    middle_alpha = top_pairing[1][0]
    middle_beta  = top_pairing[1][1]
    assert middle_alpha in ['VA','JA']
    assert middle_beta in ['VB','JB']
    other_alpha = 'JA' if middle_alpha=='VA' else 'VA'
    other_beta  = 'JB' if middle_beta =='VB' else 'VB'


    
    

    if force_pairing_order:
        assert len(force_pairing_order) == 4
        reps = force_pairing_order[:]
    else:
        reps = [ other_alpha, middle_alpha, middle_beta, other_beta ]
'''