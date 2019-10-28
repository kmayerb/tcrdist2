from tcrdist import svg_basic

def plot_pwm(SM,  my_width = 6, my_height= 6, create_file = False, output = "test.svg"):
    """
    Produces position-wise matrix (pwm) plot.

    Parameters
    ----------
    SM: StoreIOMotif
        StoreIOMotif_instance

    ts : TCRsubset

    output : str
        output svg filename

    Returns
    -------
    svg : string


    Notes
    -----

    This plot is composed of five panels.

    panel_00 : row0,col0 - V gene usage (stack of gene usage)
    panel_01 : row0,col1 - CDR3 Amino Acid Stack
    panel_02 : row0,col1 - J gene usage (stack of gene usage)
    panel_11 : row1,col1 - nucleotide source information
    panel_21 : row2,col1 - relative_entropy scaled CDR3

    The key elements pulled from the StoreIOMotif_instance in the function

        SM.ab
        SM.vl_nbr
        SM.entropy.pwm
        SM.entropy.npwm
        SM.entropy.npwm.keys
        SM.entropy.npwm
        SM.entropy.pwm, SM.entropy.scale_by_relent
        SM.jl_nbr
        SM.matched_tcrs
        SM.matched_tcrs_plus_nbrs
        SM.ep, SM.ab, len(SM.nseqs),
        SM.chi_squared
        SM.showmotif
        SM.nseqs
        SM.nseqs
        SM.expect_random float(SM.expect_nextgen)
        SM.expect_random
        SM.expect_nextgen


    """
    chain = SM.ab
    junction_bars = True
    junction_bars_order = { 'B': ['V','N1','D','N2','J'], 'A': ['V','N','J'] }
    font_family="Droid Sans Mono"
    junction_bars_color = { 'V':  'silver',
                            'N1': 'red',
                            'N':  'red',
                            'N2': 'red',
                            'D':  'black',
                            'J':  'dimgray' }
    ypad = 30
    gene_stack_width = 100
    pwm_stack_width = 300
    pwm_stack_height = 80
    svg_width = 0
    xmargin = 0
    sep = 5
    fontsize = pwm_stack_height/10
    junction_bars_ypad = 2

    cmds = []

    total_y = 0
    x0 = xmargin;
    x1 = x0 + gene_stack_width

    y0=total_y + ypad
    x0=xmargin;
    x1 = x0+gene_stack_width

    y0=total_y
    y1=y0 ## pwm
    y2=y1 ## VNDNJ pwm
    y3=y2+    pwm_stack_height ## nbr-pwm
    y4=y3+0.5*pwm_stack_height ## nbr-NVDNJ pwm
    y5=y4 ## lenseq pwm
    y6=y5 ## fwd pwm
    y7=y6 ## ng-fwd pwm
    y8=y7 ## rev pwm
    y9=y8 ## ng-rev pwm
    y10=y9+   pwm_stack_height ## nbr-pwm scaled by relent
    y_last = y10



    p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
    cmds.append( svg_basic.make_stack( p0, p1, SM.vl_nbr ) )
    cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

    x0 = x1+sep; x1 = x0+pwm_stack_width

    ## nbr pwm
    p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
    cmds.append( svg_basic.protein_logo( p0, p1, SM.entropy.pwm ) )
    cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

    if junction_bars:
        junction_pwm = SM.entropy.npwm
        junction_bars_height = y4-y3 - 2*junction_bars_ypad
        ncols = len( SM.entropy.npwm.keys() )
        junction_bar_width = ( x1-x0 )/float(ncols)
        for j in range(ncols):
            lcol = [ ( junction_pwm[j][x],x) for x in junction_bars_order[chain] ]
            y1shift = y3+junction_bars_ypad
            ## largest at the top
            for frac,a in lcol:
                y1shift_next = y1shift + frac * junction_bars_height
                color = junction_bars_color[ a ]
                p0 = [ x0+ j   *junction_bar_width, y1shift]
                p1 = [ x0+(j+1)*junction_bar_width, y1shift_next ]
                cmds.append( svg_basic.rectangle( p0, p1, fill=color, stroke=color ) )
                y1shift = y1shift_next
            assert abs( y1shift+junction_bars_ypad-y4)<1e-3
    else:
        p0,p1 = [ ( x0, y3 ), ( x1, y4 ) ]
        cmds.append( svg_basic.protein_logo( p0, p1, SM.entropy.npwm ) )
        cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

    ## nbr pwm SCALED by relent
    p0,p1 = [ ( x0, y9 ), ( x1, y10 ) ]
    cmds.append( svg_basic.protein_logo( p0, p1, SM.entropy.pwm, SM.entropy.scale_by_relent ) )
    ## box around all 3 of them
    p0,p1 = [ ( x0, y2 ), ( x1, y10 ) ]
    cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

    # RIGHT_MOST PANEL
    x0 = x1+sep; x1 = x0+gene_stack_width
    p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
    cmds.append( svg_basic.make_stack( p0, p1, SM.jl_nbr ) )
    cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

    expect_very_small = 0.001
    n= len(SM.matched_tcrs)
    nplus= len(SM.matched_tcrs_plus_nbrs)
    messages = [ '{} {} #clones={}'.format(SM.ep,SM.ab,len(SM.nseqs)),
         'chi-sq: {:.1f}'.format(float(SM.chi_squared)),
         'motif: {}'.format(''.join(SM.showmotif)),
         'match-: {} {:.1f}%'.format(n,(100.0*n)/len(SM.nseqs)),
         'match+: {} {:.1f}%'.format(nplus,(100.0*nplus)/len(SM.nseqs)),
         'expect: {:.3f}'.format(max(float(SM.expect_random),float(SM.expect_nextgen))),
         'enrich: {:.1f}'.format(float(n)/max([float(SM.expect_random),
                                               float(SM.expect_nextgen),
                                               float(expect_very_small)])) ]

    for ii,text in enumerate(messages):
        p0 = ( x0, y4 + (ii+1)*fontsize ) ## lower left of text

        width = fontsize * len(text) * 0.6
        svg_width = max( svg_width, p0[0] + width + 50 )

        cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )
    cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

    y_last = y10
    if create_file:
        svg_basic.create_file( cmds, x1, y_last, output, create_png=False )
    else:
        svg_txt = create_svg(cmds, width = my_width, height= my_height)
        return(svg_txt)

def create_svg( cmds, width, height, background_color=None, use_xlink=False ):
    out = ''
    extra = '' if not use_xlink else 'xmlns:xlink="http://www.w3.org/1999/xlink"'
    out += '<svg width="{}" height="{}" xmlns="http://www.w3.org/2000/svg" version="1.1" {} >\n'\
              .format(int(width),int(height),extra)
    if background_color:
        out += rectangle( (0,0), (width, height), background_color, 'white', 0 )
    out += '\n'.join(cmds)+'\n'
    out += '</svg>\n'
    return out
