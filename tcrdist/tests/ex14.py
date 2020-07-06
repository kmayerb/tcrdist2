import pytest 
import tcrsampler 


def test_ex14():
    import pandas as pd
    import tcrsampler as ts
    from tcrdist.repertoire import TCRrep
    from tcrdist import cluster
    import pwseqdist as pw


    df = pd.read_csv('dash.csv')
    df = df[df.epitope.isin(['PA'])]
    tr = TCRrep(cell_df=df, chains=['alpha','beta'], organism='mouse')
    tr.tcrdist2(processes = 6,
                metric = 'nw',
                reduce = True,
                dump = False,
                save = False)

    import tcrsampler as ts
    
    # select a mouse and beta chain specific background from ruggiero et al. 2015
    tsbeta = ts.TCRsampler(default_background = 'ruggiero_mouse_beta_t.tsv.sampler.tsv')
    # create a deeper background set. 
    tsbeta.build_background(max_rows=1000)
    
    #t = TCRsampler(use_default=True)

    tr.cluster_index = tr.simple_cluster_index(pw_distances = tr.pw_beta,
                                               method = 'complete',
                                               criterion = "distance",
                                               t = 50)

    assert len(tr.cluster_index) == tr.clone_df.shape[0]
    tr.cluster_df = tr.cluster_index_to_df(tr.cluster_index)
    
    motifs = tr.generate_cluster_summary(sampler = tsbeta, n = 200)

    # you can visualize clusters wiht a simple index.html
    
    # <!DOCTYPE html>
    # <html>
    #   <frameset cols="50%,50%" marginwidth="20">
    #     <frame src="static/toc.html">
    #     <frame src="static/5_cluster_id.html"name="logo_frame">
    #   </frameset>
    # </html>
    
    # # all the motifs can be inspected by cluster number
    # motifs[2].cluster_attributes
    # motifs[2].pal_stat
    # motifs[2].pal_motif
    
    
    
    # all_pal_svgs = [motifs[k].pal_svg for k in sorted(tr.cluster_df.cluster_id)]
    # #cluster._write_svgs(svgs = all_pal_svgs, name = "all_pal_svgs.svg", dest = 'static')
    # cluster._write_svgs_w_labs(svgs = all_pal_svgs, labs = sorted(tr.cluster_df.cluster_id), name = "all_pal_svgs.svg", dest = 'static')






    


