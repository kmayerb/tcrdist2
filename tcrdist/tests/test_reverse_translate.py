import pytest
from tcrdist.reverse_translate import TCRcodon   
from tcrdist.tests.my_test_subset import clone_df_subset                                                                                                                                                                                                                                                                                                           
import numpy as np
import pandas as pd    
 
def test_TCRcodon_mouse_ab_init():
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    assert set(tc.all_genes.keys()) == set(['TRAV1*01', 'TRAV1*02', 'TRAV10*01', 'TRAV10*02', 'TRAV10*03', 'TRAV10*04', 'TRAV10*05', 'TRAV10D*01', 'TRAV10D*02', 'TRAV10N*01', 'TRAV11*01', 'TRAV11*02', 'TRAV11D*01', 'TRAV11N*01', 'TRAV12-1*01', 'TRAV12-1*02', 'TRAV12-1*03', 'TRAV12-1*04', 'TRAV12-1*05', 'TRAV12-2*01', 'TRAV12-3*01', 'TRAV12-3*02', 'TRAV12-3*03', 'TRAV12-3*04', 'TRAV12D-1*01', 'TRAV12D-1*02', 'TRAV12D-1*03', 'TRAV12D-1*04', 'TRAV12D-1*05', 'TRAV12D-2*01', 'TRAV12D-2*02', 'TRAV12D-2*03', 'TRAV12D-2*04', 'TRAV12D-2*05', 'TRAV12D-3*01', 'TRAV12D-3*02', 'TRAV12D-3*03', 'TRAV12N-1*01', 'TRAV12N-2*01', 'TRAV12N-3*01', 'TRAV13-1*01', 'TRAV13-2*01', 'TRAV13-2*02', 'TRAV13-3*01', 'TRAV13-3*02', 'TRAV13-4/DV7*01', 'TRAV13-4/DV7*02', 'TRAV13-4/DV7*03', 'TRAV13-5*01', 'TRAV13D-1*01', 'TRAV13D-1*02', 'TRAV13D-1*03', 'TRAV13D-2*01', 'TRAV13D-2*02', 'TRAV13D-3*01', 'TRAV13D-3*02', 'TRAV13D-4*01', 'TRAV13D-4*02', 'TRAV13D-4*03', 'TRAV13N-1*01', 'TRAV13N-2*01', 'TRAV13N-3*01', 'TRAV13N-4*01', 'TRAV14-1*01', 'TRAV14-1*02', 'TRAV14-1*03', 'TRAV14-2*01', 'TRAV14-2*02', 'TRAV14-2*03', 'TRAV14-3*01', 'TRAV14-3*02', 'TRAV14-3*03', 'TRAV14D-1*01', 'TRAV14D-1*02', 'TRAV14D-2*01', 'TRAV14D-2*02', 'TRAV14D-2*03', 'TRAV14D-3/DV8*01', 'TRAV14D-3/DV8*02', 'TRAV14D-3/DV8*03', 'TRAV14D-3/DV8*04', 'TRAV14D-3/DV8*05', 'TRAV14D-3/DV8*06', 'TRAV14D-3/DV8*07', 'TRAV14D-3/DV8*08', 'TRAV14N-1*01', 'TRAV14N-2*01', 'TRAV14N-3*01', 'TRAV15-1/DV6-1*01', 'TRAV15-1/DV6-1*02', 'TRAV15-2/DV6-2*01', 'TRAV15-2/DV6-2*02', 'TRAV15D-1/DV6D-1*01', 'TRAV15D-1/DV6D-1*02', 'TRAV15D-1/DV6D-1*03', 'TRAV15D-1/DV6D-1*04', 'TRAV15D-1/DV6D-1*05', 'TRAV15D-1/DV6D-1*06', 'TRAV15D-2/DV6D-2*01', 'TRAV15D-2/DV6D-2*02', 'TRAV15D-2/DV6D-2*03', 'TRAV15D-2/DV6D-2*04', 'TRAV15D-2/DV6D-2*05', 'TRAV15N-1*01', 'TRAV15N-2*01', 'TRAV16*01', 'TRAV16*02', 'TRAV16*03', 'TRAV16*04', 'TRAV16*05', 'TRAV16D/DV11*01', 'TRAV16D/DV11*02', 'TRAV16D/DV11*03', 'TRAV16N*01', 'TRAV17*01', 'TRAV17*02', 'TRAV18*01', 'TRAV19*01', 'TRAV19*03', 'TRAV2*01', 'TRAV20*01', 'TRAV20*02', 'TRAV21/DV12*01', 'TRAV21/DV12*02', 'TRAV3-1*01', 'TRAV3-1*02', 'TRAV3-3*01', 'TRAV3-4*01', 'TRAV3D-3*01', 'TRAV3D-3*02', 'TRAV3N-3*01', 'TRAV4-2*01', 'TRAV4-2*02', 'TRAV4-3*01', 'TRAV4-3*02', 'TRAV4-4/DV10*01', 'TRAV4D-2*01', 'TRAV4D-3*01', 'TRAV4D-3*02', 'TRAV4D-3*03', 'TRAV4D-3*04', 'TRAV4D-4*01', 'TRAV4D-4*02', 'TRAV4D-4*03', 'TRAV4D-4*04', 'TRAV4N-3*01', 'TRAV4N-4*01', 'TRAV5-1*01', 'TRAV5-2*01', 'TRAV5-4*01', 'TRAV5D-2*01', 'TRAV5D-4*01', 'TRAV5D-4*02', 'TRAV5D-4*03', 'TRAV5D-4*04', 'TRAV5D-4*05', 'TRAV5N-2*01', 'TRAV5N-4*01', 'TRAV6-1*01', 'TRAV6-1*02', 'TRAV6-2*01', 'TRAV6-2*02', 'TRAV6-2*03', 'TRAV6-3*01', 'TRAV6-3*02', 'TRAV6-4*01', 'TRAV6-4*02', 'TRAV6-4*03', 'TRAV6-5*01', 'TRAV6-5*02', 'TRAV6-5*03', 'TRAV6-5*04', 'TRAV6-6*01', 'TRAV6-6*02', 'TRAV6-6*03', 'TRAV6-7/DV9*01', 'TRAV6-7/DV9*02', 'TRAV6-7/DV9*03', 'TRAV6-7/DV9*04', 'TRAV6-7/DV9*06', 'TRAV6-7/DV9*07', 'TRAV6-7/DV9*08', 'TRAV6D-3*01', 'TRAV6D-3*02', 'TRAV6D-4*01', 'TRAV6D-5*01', 'TRAV6D-6*01', 'TRAV6D-6*02', 'TRAV6D-6*03', 'TRAV6D-6*04', 'TRAV6D-6*05', 'TRAV6D-7*01', 'TRAV6D-7*02', 'TRAV6D-7*03', 'TRAV6D-7*04', 'TRAV6N-5*01', 'TRAV6N-6*01', 'TRAV6N-7*01', 'TRAV7-1*01', 'TRAV7-2*01', 'TRAV7-2*02', 'TRAV7-3*01', 'TRAV7-3*02', 'TRAV7-3*03', 'TRAV7-3*04', 'TRAV7-4*01', 'TRAV7-4*02', 'TRAV7-5*01', 'TRAV7-5*02', 'TRAV7-5*03', 'TRAV7-6*01', 'TRAV7-6*02', 'TRAV7D-2*01', 'TRAV7D-2*02', 'TRAV7D-2*03', 'TRAV7D-3*01', 'TRAV7D-3*02', 'TRAV7D-4*01', 'TRAV7D-4*02', 'TRAV7D-4*03', 'TRAV7D-5*01', 'TRAV7D-6*01', 'TRAV7D-6*02', 'TRAV7N-4*01', 'TRAV7N-5*01', 'TRAV7N-6*01', 'TRAV8-1*01', 'TRAV8-1*02', 'TRAV8-1*03', 'TRAV8-2*01', 'TRAV8D-1*01', 'TRAV8D-1*02', 'TRAV8D-2*01', 'TRAV8D-2*02', 'TRAV8D-2*03', 'TRAV8N-2*01', 'TRAV9-1*01', 'TRAV9-1*02', 'TRAV9-2*01', 'TRAV9-3*01', 'TRAV9-3*02', 'TRAV9-3*03', 'TRAV9-4*01', 'TRAV9D-1*01', 'TRAV9D-1*02', 'TRAV9D-2*01', 'TRAV9D-2*02', 'TRAV9D-2*03', 'TRAV9D-3*01', 'TRAV9D-3*02', 'TRAV9D-4*01', 'TRAV9D-4*03', 'TRAV9D-4*04', 'TRAV9N-2*01', 'TRAV9N-3*01', 'TRAV9N-4*01', 'TRAJ11*01', 'TRAJ12*01', 'TRAJ13*01', 'TRAJ15*01', 'TRAJ16*01', 'TRAJ17*01', 'TRAJ18*01', 'TRAJ19*01', 'TRAJ2*01', 'TRAJ2*02', 'TRAJ20*01', 'TRAJ21*01', 'TRAJ21*02', 'TRAJ22*01', 'TRAJ23*01', 'TRAJ24*01', 'TRAJ25*01', 'TRAJ26*01', 'TRAJ27*01', 'TRAJ28*01', 'TRAJ29*01', 'TRAJ3*01', 'TRAJ30*01', 'TRAJ31*01', 'TRAJ31*02', 'TRAJ32*01', 'TRAJ32*02', 'TRAJ33*01', 'TRAJ34*01', 'TRAJ34*02', 'TRAJ35*01', 'TRAJ35*02', 'TRAJ36*01', 'TRAJ37*01', 'TRAJ38*01', 'TRAJ39*01', 'TRAJ4*01', 'TRAJ40*01', 'TRAJ41*01', 'TRAJ42*01', 'TRAJ42*02', 'TRAJ43*01', 'TRAJ43*02', 'TRAJ44*01', 'TRAJ45*01', 'TRAJ45*02', 'TRAJ46*01', 'TRAJ47*01', 'TRAJ48*01', 'TRAJ49*01', 'TRAJ5*01', 'TRAJ50*01', 'TRAJ52*01', 'TRAJ53*01', 'TRAJ54*01', 'TRAJ56*01', 'TRAJ57*01', 'TRAJ58*01', 'TRAJ59*01', 'TRAJ6*01', 'TRAJ60*01', 'TRAJ61*01', 'TRAJ7*01', 'TRAJ9*01', 'TRAJ9*02', 'TRBV1*01', 'TRBV1*02', 'TRBV10*01', 'TRBV12-1*01', 'TRBV12-1*02', 'TRBV12-2*01', 'TRBV12-2*02', 'TRBV13-1*01', 'TRBV13-1*02', 'TRBV13-2*01', 'TRBV13-2*02', 'TRBV13-2*03', 'TRBV13-2*04', 'TRBV13-2*05', 'TRBV13-3*01', 'TRBV14*01', 'TRBV15*01', 'TRBV16*01', 'TRBV16*02', 'TRBV16*03', 'TRBV16*04', 'TRBV17*01', 'TRBV19*01', 'TRBV19*02', 'TRBV19*03', 'TRBV2*01', 'TRBV20*01', 'TRBV20*02', 'TRBV21*01', 'TRBV23*01', 'TRBV24*01', 'TRBV24*02', 'TRBV24*03', 'TRBV24*04', 'TRBV26*01', 'TRBV26*02', 'TRBV29*01', 'TRBV29*02', 'TRBV3*01', 'TRBV3*02', 'TRBV30*01', 'TRBV31*01', 'TRBV31*02', 'TRBV4*01', 'TRBV4*02', 'TRBV5*01', 'TRBV5*02', 'TRBV5*03', 'TRBV5*04', 'TRBV5*05', 'TRBV8*01', 'TRBV9*01', 'TRBJ1-1*01', 'TRBJ1-1*02', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-4*02', 'TRBJ1-5*01', 'TRBJ1-5*02', 'TRBJ1-5*03', 'TRBJ1-6*01', 'TRBJ1-7*01', 'TRBJ2-1*01', 'TRBJ2-2*01', 'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-6*01', 'TRBJ2-7*01', 'TRBJ2-7*02', 'TRBD1*01', 'TRBD2*01'])

def test_TCRcodon_human_ab_init():
    tc = TCRcodon(organism = "human", db_file = "alphabeta_db.tsv")
    assert set(tc.all_genes.keys()) ==  set(['TRAV1-1*01', 'TRAV1-1*02', 'TRAV1-2*01', 'TRAV1-2*02', 'TRAV10*01', 'TRAV11*01', 'TRAV12-1*01', 'TRAV12-1*02', 'TRAV12-2*01', 'TRAV12-2*02', 'TRAV12-2*03', 'TRAV12-3*01', 'TRAV12-3*02', 'TRAV13-1*01', 'TRAV13-1*02', 'TRAV13-1*03', 'TRAV13-2*01', 'TRAV13-2*02', 'TRAV14/DV4*01', 'TRAV14/DV4*02', 'TRAV14/DV4*03', 'TRAV14/DV4*04', 'TRAV16*01', 'TRAV17*01', 'TRAV18*01', 'TRAV19*01', 'TRAV2*01', 'TRAV2*02', 'TRAV20*01', 'TRAV20*02', 'TRAV20*03', 'TRAV20*04', 'TRAV21*01', 'TRAV21*02', 'TRAV22*01', 'TRAV23/DV6*01', 'TRAV23/DV6*02', 'TRAV23/DV6*03', 'TRAV23/DV6*04', 'TRAV24*01', 'TRAV24*02', 'TRAV25*01', 'TRAV26-1*01', 'TRAV26-1*02', 'TRAV26-1*03', 'TRAV26-2*01', 'TRAV26-2*02', 'TRAV27*01', 'TRAV27*02', 'TRAV27*03', 'TRAV29/DV5*01', 'TRAV29/DV5*02', 'TRAV3*01', 'TRAV30*01', 'TRAV30*02', 'TRAV30*03', 'TRAV30*04', 'TRAV34*01', 'TRAV35*01', 'TRAV35*02', 'TRAV36/DV7*01', 'TRAV36/DV7*02', 'TRAV36/DV7*03', 'TRAV36/DV7*04', 'TRAV38-1*01', 'TRAV38-1*02', 'TRAV38-1*03', 'TRAV38-1*04', 'TRAV38-2/DV8*01', 'TRAV39*01', 'TRAV4*01', 'TRAV40*01', 'TRAV41*01', 'TRAV5*01', 'TRAV6*01', 'TRAV6*02', 'TRAV6*03', 'TRAV6*04', 'TRAV6*05', 'TRAV6*06', 'TRAV7*01', 'TRAV8-1*01', 'TRAV8-1*02', 'TRAV8-2*01', 'TRAV8-2*02', 'TRAV8-3*01', 'TRAV8-3*02', 'TRAV8-3*03', 'TRAV8-4*01', 'TRAV8-4*02', 'TRAV8-4*03', 'TRAV8-4*04', 'TRAV8-4*05', 'TRAV8-4*06', 'TRAV8-4*07', 'TRAV8-6*01', 'TRAV8-6*02', 'TRAV8-7*01', 'TRAV9-1*01', 'TRAV9-2*01', 'TRAV9-2*02', 'TRAV9-2*03', 'TRAV9-2*04', 'TRAJ1*01', 'TRAJ10*01', 'TRAJ11*01', 'TRAJ12*01', 'TRAJ13*01', 'TRAJ13*02', 'TRAJ14*01', 'TRAJ15*01', 'TRAJ15*02', 'TRAJ16*01', 'TRAJ17*01', 'TRAJ18*01', 'TRAJ19*01', 'TRAJ2*01', 'TRAJ20*01', 'TRAJ21*01', 'TRAJ22*01', 'TRAJ23*01', 'TRAJ23*02', 'TRAJ24*01', 'TRAJ24*02', 'TRAJ25*01', 'TRAJ26*01', 'TRAJ27*01', 'TRAJ28*01', 'TRAJ29*01', 'TRAJ3*01', 'TRAJ30*01', 'TRAJ31*01', 'TRAJ32*01', 'TRAJ32*02', 'TRAJ33*01', 'TRAJ34*01', 'TRAJ35*01', 'TRAJ36*01', 'TRAJ37*01', 'TRAJ37*02', 'TRAJ38*01', 'TRAJ39*01', 'TRAJ4*01', 'TRAJ40*01', 'TRAJ41*01', 'TRAJ42*01', 'TRAJ43*01', 'TRAJ44*01', 'TRAJ45*01', 'TRAJ46*01', 'TRAJ47*01', 'TRAJ47*02', 'TRAJ48*01', 'TRAJ49*01', 'TRAJ5*01', 'TRAJ50*01', 'TRAJ51*01', 'TRAJ52*01', 'TRAJ53*01', 'TRAJ54*01', 'TRAJ55*01', 'TRAJ56*01', 'TRAJ57*01', 'TRAJ58*01', 'TRAJ59*01', 'TRAJ6*01', 'TRAJ60*01', 'TRAJ61*01', 'TRAJ7*01', 'TRAJ8*01', 'TRAJ9*01', 'TRBV1*01', 'TRBV10-1*01', 'TRBV10-1*02', 'TRBV10-2*01', 'TRBV10-2*02', 'TRBV10-3*01', 'TRBV10-3*02', 'TRBV10-3*03', 'TRBV10-3*04', 'TRBV11-1*01', 'TRBV11-2*01', 'TRBV11-2*02', 'TRBV11-2*03', 'TRBV11-3*01', 'TRBV11-3*02', 'TRBV11-3*03', 'TRBV12-1*01', 'TRBV12-2*01', 'TRBV12-3*01', 'TRBV12-4*01', 'TRBV12-4*02', 'TRBV12-5*01', 'TRBV13*01', 'TRBV13*02', 'TRBV14*01', 'TRBV14*02', 'TRBV15*01', 'TRBV15*02', 'TRBV15*03', 'TRBV16*01', 'TRBV16*02', 'TRBV16*03', 'TRBV17*01', 'TRBV18*01', 'TRBV19*01', 'TRBV19*02', 'TRBV19*03', 'TRBV2*01', 'TRBV2*02', 'TRBV2*03', 'TRBV20-1*01', 'TRBV20-1*02', 'TRBV20-1*03', 'TRBV20-1*04', 'TRBV20-1*05', 'TRBV20-1*06', 'TRBV20-1*07', 'TRBV20/OR9-2*01', 'TRBV20/OR9-2*02', 'TRBV20/OR9-2*03', 'TRBV21-1*01', 'TRBV21/OR9-2*01', 'TRBV23-1*01', 'TRBV23/OR9-2*01', 'TRBV23/OR9-2*02', 'TRBV24-1*01', 'TRBV24/OR9-2*01', 'TRBV25-1*01', 'TRBV26*01', 'TRBV26/OR9-2*01', 'TRBV26/OR9-2*02', 'TRBV27*01', 'TRBV28*01', 'TRBV29-1*01', 'TRBV29-1*02', 'TRBV29-1*03', 'TRBV29/OR9-2*01', 'TRBV29/OR9-2*02', 'TRBV3-1*01', 'TRBV3-1*02', 'TRBV3-2*01', 'TRBV3-2*02', 'TRBV3-2*03', 'TRBV30*01', 'TRBV30*02', 'TRBV30*04', 'TRBV30*05', 'TRBV4-1*01', 'TRBV4-1*02', 'TRBV4-2*01', 'TRBV4-2*02', 'TRBV4-3*01', 'TRBV4-3*02', 'TRBV4-3*03', 'TRBV4-3*04', 'TRBV5-1*01', 'TRBV5-1*02', 'TRBV5-3*01', 'TRBV5-3*02', 'TRBV5-4*01', 'TRBV5-4*02', 'TRBV5-4*03', 'TRBV5-4*04', 'TRBV5-5*01', 'TRBV5-5*02', 'TRBV5-5*03', 'TRBV5-6*01', 'TRBV5-7*01', 'TRBV5-8*01', 'TRBV5-8*02', 'TRBV6-1*01', 'TRBV6-2*01', 'TRBV6-3*01', 'TRBV6-4*01', 'TRBV6-4*02', 'TRBV6-5*01', 'TRBV6-6*01', 'TRBV6-6*02', 'TRBV6-6*03', 'TRBV6-6*04', 'TRBV6-6*05', 'TRBV6-7*01', 'TRBV6-8*01', 'TRBV6-9*01', 'TRBV7-1*01', 'TRBV7-2*01', 'TRBV7-2*02', 'TRBV7-2*03', 'TRBV7-2*04', 'TRBV7-3*01', 'TRBV7-3*02', 'TRBV7-3*03', 'TRBV7-3*04', 'TRBV7-3*05', 'TRBV7-4*01', 'TRBV7-6*01', 'TRBV7-6*02', 'TRBV7-7*01', 'TRBV7-7*02', 'TRBV7-8*01', 'TRBV7-8*02', 'TRBV7-8*03', 'TRBV7-9*01', 'TRBV7-9*02', 'TRBV7-9*03', 'TRBV7-9*04', 'TRBV7-9*05', 'TRBV7-9*06', 'TRBV7-9*07', 'TRBV9*01', 'TRBV9*02', 'TRBV9*03', 'TRBJ1-1*01', 'TRBJ1-2*01', 'TRBJ1-3*01', 'TRBJ1-4*01', 'TRBJ1-5*01', 'TRBJ1-6*01', 'TRBJ1-6*02', 'TRBJ2-1*01', 'TRBJ2-2*01', 'TRBJ2-2P*01', 'TRBJ2-3*01', 'TRBJ2-4*01', 'TRBJ2-5*01', 'TRBJ2-6*01', 'TRBJ2-7*01', 'TRBJ2-7*02', 'TRBD1*01', 'TRBD2*01', 'TRBD2*02'])

def test_TCRcodon_mouse_gd_init():
    tc = TCRcodon(organism = "mouse", db_file = "gammadelta_db.tsv")
    assert set(tc.all_genes.keys()) == set(['TRGV1*01', 'TRGV1*02', 'TRGV1*03', 'TRGV1*04', 'TRGV1*05', 'TRGV1*06', 'TRGV1*07', 'TRGV1*08', 'TRGV2*01', 'TRGV2*02', 'TRGV2*03', 'TRGV2*04', 'TRGV2*05', 'TRGV3*01', 'TRGV3*02', 'TRGV3*03', 'TRGV4*01', 'TRGV4*02', 'TRGV4*03', 'TRGV4*04', 'TRGV4*05', 'TRGV5*01', 'TRGV6*01', 'TRGV6*02', 'TRGV6*03', 'TRGV6*04', 'TRGV7*01', 'TRGV7*02', 'TRGJ3*01', 'TRGJ2*01', 'TRGJ1*01', 'TRGJ4*01', 'TRAV1*01', 'TRAV1*02', 'TRAV10*01', 'TRAV10*02', 'TRAV10*03', 'TRAV10*04', 'TRAV10*05', 'TRAV10D*01', 'TRAV10D*02', 'TRAV10N*01', 'TRAV11*01', 'TRAV11*02', 'TRAV11D*01', 'TRAV11N*01', 'TRAV12-1*01', 'TRAV12-1*02', 'TRAV12-1*03', 'TRAV12-1*04', 'TRAV12-1*05', 'TRAV12-2*01', 'TRAV12-3*01', 'TRAV12-3*02', 'TRAV12-3*03', 'TRAV12-3*04', 'TRAV12D-1*01', 'TRAV12D-1*02', 'TRAV12D-1*04', 'TRAV12D-1*05', 'TRAV12D-2*01', 'TRAV12D-2*02', 'TRAV12D-2*03', 'TRAV12D-2*04', 'TRAV12D-2*05', 'TRAV12D-3*01', 'TRAV12D-3*02', 'TRAV12D-3*03', 'TRAV12N-1*01', 'TRAV12N-2*01', 'TRAV12N-3*01', 'TRAV13-1*01', 'TRAV13-2*01', 'TRAV13-2*02', 'TRAV13-3*01', 'TRAV13-3*02', 'TRAV13-4/DV7*01', 'TRAV13-4/DV7*02', 'TRAV13-4/DV7*03', 'TRAV13-5*01', 'TRAV13D-1*01', 'TRAV13D-1*02', 'TRAV13D-1*03', 'TRAV13D-2*01', 'TRAV13D-2*02', 'TRAV13D-3*01', 'TRAV13D-3*02', 'TRAV13D-4*01', 'TRAV13D-4*02', 'TRAV13D-4*03', 'TRAV13N-1*01', 'TRAV13N-2*01', 'TRAV13N-3*01', 'TRAV13N-4*01', 'TRAV14-1*01', 'TRAV14-1*02', 'TRAV14-1*03', 'TRAV14-2*01', 'TRAV14-2*02', 'TRAV14-2*03', 'TRAV14-3*01', 'TRAV14-3*02', 'TRAV14-3*03', 'TRAV14D-1*01', 'TRAV14D-1*02', 'TRAV14D-2*01', 'TRAV14D-2*02', 'TRAV14D-2*03', 'TRAV14D-3/DV8*01', 'TRAV14D-3/DV8*02', 'TRAV14D-3/DV8*03', 'TRAV14D-3/DV8*04', 'TRAV14D-3/DV8*05', 'TRAV14D-3/DV8*06', 'TRAV14D-3/DV8*07', 'TRAV14D-3/DV8*08', 'TRAV14N-1*01', 'TRAV14N-2*01', 'TRAV14N-3*01', 'TRAV15-1/DV6-1*01', 'TRAV15-1/DV6-1*02', 'TRAV15-2/DV6-2*01', 'TRAV15-2/DV6-2*02', 'TRAV15D-1/DV6D-1*01', 'TRAV15D-1/DV6D-1*02', 'TRAV15D-1/DV6D-1*03', 'TRAV15D-1/DV6D-1*04', 'TRAV15D-1/DV6D-1*05', 'TRAV15D-1/DV6D-1*06', 'TRAV15D-2/DV6D-2*01', 'TRAV15D-2/DV6D-2*02', 'TRAV15D-2/DV6D-2*03', 'TRAV15D-2/DV6D-2*04', 'TRAV15D-2/DV6D-2*05', 'TRAV15N-1*01', 'TRAV15N-2*01', 'TRAV16*01', 'TRAV16*02', 'TRAV16*03', 'TRAV16*04', 'TRAV16*05', 'TRAV16D/DV11*01', 'TRAV16D/DV11*02', 'TRAV16D/DV11*03', 'TRAV16N*01', 'TRAV17*01', 'TRAV17*02', 'TRAV18*01', 'TRAV19*01', 'TRAV19*03', 'TRAV2*01', 'TRAV20*01', 'TRAV20*02', 'TRAV21/DV12*01', 'TRAV21/DV12*02', 'TRAV3-1*01', 'TRAV3-1*02', 'TRAV3-3*01', 'TRAV3-4*01', 'TRAV3D-3*01', 'TRAV3D-3*02', 'TRAV3N-3*01', 'TRAV4-2*01', 'TRAV4-2*02', 'TRAV4-3*01', 'TRAV4-3*02', 'TRAV4-4/DV10*01', 'TRAV4D-2*01', 'TRAV4D-3*01', 'TRAV4D-3*02', 'TRAV4D-3*03', 'TRAV4D-3*04', 'TRAV4D-4*01', 'TRAV4D-4*02', 'TRAV4D-4*03', 'TRAV4D-4*04', 'TRAV4N-3*01', 'TRAV4N-4*01', 'TRAV5-1*01', 'TRAV5-2*01', 'TRAV5-4*01', 'TRAV5D-2*01', 'TRAV5D-4*01', 'TRAV5D-4*02', 'TRAV5D-4*03', 'TRAV5D-4*04', 'TRAV5D-4*05', 'TRAV5N-2*01', 'TRAV5N-4*01', 'TRAV6-1*01', 'TRAV6-1*02', 'TRAV6-2*01', 'TRAV6-2*02', 'TRAV6-2*03', 'TRAV6-3*01', 'TRAV6-3*02', 'TRAV6-4*01', 'TRAV6-4*02', 'TRAV6-4*03', 'TRAV6-5*01', 'TRAV6-5*02', 'TRAV6-5*03', 'TRAV6-5*04', 'TRAV6-6*01', 'TRAV6-6*02', 'TRAV6-6*03', 'TRAV6-7/DV9*01', 'TRAV6-7/DV9*02', 'TRAV6-7/DV9*03', 'TRAV6-7/DV9*04', 'TRAV6-7/DV9*06', 'TRAV6-7/DV9*07', 'TRAV6-7/DV9*08', 'TRAV6D-3*01', 'TRAV6D-3*02', 'TRAV6D-4*01', 'TRAV6D-5*01', 'TRAV6D-6*01', 'TRAV6D-6*02', 'TRAV6D-6*03', 'TRAV6D-6*04', 'TRAV6D-6*05', 'TRAV6D-7*01', 'TRAV6D-7*02', 'TRAV6D-7*03', 'TRAV6D-7*04', 'TRAV6N-5*01', 'TRAV6N-6*01', 'TRAV6N-7*01', 'TRAV7-1*01', 'TRAV7-2*01', 'TRAV7-2*02', 'TRAV7-3*01', 'TRAV7-3*02', 'TRAV7-3*03', 'TRAV7-3*04', 'TRAV7-4*01', 'TRAV7-4*02', 'TRAV7-5*01', 'TRAV7-5*02', 'TRAV7-5*03', 'TRAV7-6*01', 'TRAV7-6*02', 'TRAV7D-2*01', 'TRAV7D-2*02', 'TRAV7D-2*03', 'TRAV7D-3*01', 'TRAV7D-3*02', 'TRAV7D-4*01', 'TRAV7D-4*02', 'TRAV7D-4*03', 'TRAV7D-5*01', 'TRAV7D-6*01', 'TRAV7D-6*02', 'TRAV7N-4*01', 'TRAV7N-5*01', 'TRAV7N-6*01', 'TRAV8-1*01', 'TRAV8-1*02', 'TRAV8-1*03', 'TRAV8-2*01', 'TRAV8D-1*01', 'TRAV8D-1*02', 'TRAV8D-2*01', 'TRAV8D-2*02', 'TRAV8D-2*03', 'TRAV8N-2*01', 'TRAV9-1*01', 'TRAV9-1*02', 'TRAV9-2*01', 'TRAV9-3*01', 'TRAV9-3*02', 'TRAV9-3*03', 'TRAV9-4*01', 'TRAV9D-1*01', 'TRAV9D-1*02', 'TRAV9D-2*01', 'TRAV9D-2*02', 'TRAV9D-2*03', 'TRAV9D-3*01', 'TRAV9D-3*02', 'TRAV9D-4*01', 'TRAV9D-4*03', 'TRAV9D-4*04', 'TRAV9N-2*01', 'TRAV9N-3*01', 'TRAV9N-4*01', 'TRDV1*01', 'TRDV2-1*01', 'TRDV2-2*01', 'TRDV2-2*02', 'TRDV4*01', 'TRDV5*01', 'TRDV5*02', 'TRDV5*03', 'TRDV5*04', 'TRDJ1*01', 'TRDJ2*01', 'TRDJ2*02', 'TRDD1*01', 'TRDD2*01'])

def test_TCRcodon_human_gd_init():
    tc = TCRcodon(organism = "human", db_file = "gammadelta_db.tsv")
    assert set(tc.all_genes.keys()) == set(['TRGV1*01', 'TRGV10*01', 'TRGV10*02', 'TRGV11*01', 'TRGV11*02', 'TRGV2*01', 'TRGV2*02', 'TRGV3*01', 'TRGV3*02', 'TRGV4*01', 'TRGV4*02', 'TRGV5*01', 'TRGV5P*01', 'TRGV5P*02', 'TRGV8*01', 'TRGV9*01', 'TRGV9*02', 'TRGVA*01', 'TRGJ1*01', 'TRGJ1*02', 'TRGJP1*01', 'TRGJ2*01', 'TRGJP*01', 'TRGJP2*01', 'TRAV1-1*01', 'TRAV1-1*02', 'TRAV1-2*01', 'TRAV1-2*02', 'TRAV10*01', 'TRAV11*01', 'TRAV12-1*01', 'TRAV12-1*02', 'TRAV12-2*01', 'TRAV12-2*02', 'TRAV12-2*03', 'TRAV12-3*01', 'TRAV12-3*02', 'TRAV13-1*01', 'TRAV13-1*02', 'TRAV13-1*03', 'TRAV13-2*01', 'TRAV13-2*02', 'TRAV14/DV4*01', 'TRAV14/DV4*02', 'TRAV14/DV4*03', 'TRAV14/DV4*04', 'TRAV16*01', 'TRAV17*01', 'TRAV18*01', 'TRAV19*01', 'TRAV2*01', 'TRAV2*02', 'TRAV20*01', 'TRAV20*02', 'TRAV20*03', 'TRAV20*04', 'TRAV21*01', 'TRAV21*02', 'TRAV22*01', 'TRAV23/DV6*01', 'TRAV23/DV6*02', 'TRAV23/DV6*03', 'TRAV23/DV6*04', 'TRAV24*01', 'TRAV24*02', 'TRAV25*01', 'TRAV26-1*01', 'TRAV26-1*02', 'TRAV26-1*03', 'TRAV26-2*01', 'TRAV26-2*02', 'TRAV27*01', 'TRAV27*02', 'TRAV27*03', 'TRAV29/DV5*01', 'TRAV29/DV5*02', 'TRAV3*01', 'TRAV30*01', 'TRAV30*02', 'TRAV30*03', 'TRAV30*04', 'TRAV34*01', 'TRAV35*01', 'TRAV35*02', 'TRAV36/DV7*01', 'TRAV36/DV7*02', 'TRAV36/DV7*03', 'TRAV36/DV7*04', 'TRAV38-1*01', 'TRAV38-1*02', 'TRAV38-1*03', 'TRAV38-1*04', 'TRAV38-2/DV8*01', 'TRAV39*01', 'TRAV4*01', 'TRAV40*01', 'TRAV41*01', 'TRAV5*01', 'TRAV6*01', 'TRAV6*02', 'TRAV6*03', 'TRAV6*04', 'TRAV6*05', 'TRAV6*06', 'TRAV7*01', 'TRAV8-1*01', 'TRAV8-1*02', 'TRAV8-2*01', 'TRAV8-2*02', 'TRAV8-3*01', 'TRAV8-3*02', 'TRAV8-3*03', 'TRAV8-4*01', 'TRAV8-4*02', 'TRAV8-4*03', 'TRAV8-4*04', 'TRAV8-4*05', 'TRAV8-4*06', 'TRAV8-4*07', 'TRAV8-6*01', 'TRAV8-6*02', 'TRAV8-7*01', 'TRAV9-1*01', 'TRAV9-2*01', 'TRAV9-2*02', 'TRAV9-2*03', 'TRAV9-2*04', 'TRDV1*01', 'TRDV2*01', 'TRDV2*02', 'TRDV2*03', 'TRDV3*01', 'TRDV3*02', 'TRDJ1*01', 'TRDJ4*01', 'TRDJ3*01', 'TRDJ2*01', 'TRDD1*01', 'TRDD3*01', 'TRDD2*01'])

def test_TCRcodon_single_example():
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    r = tc.guess_reverse_translation(v_gene_name= 'TRBV29*01' , j_gene_name= 'TRBJ1-5*01' , cdr3_aa = 'CASSEGEAPLF')
    assert r == 'TGTGCTAGCAGTGAGGGAGAGGCTCCGCTTTTT'

def test_TCRcodon_single_example2():
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    r = tc.guess_reverse_translation(v_gene_name= 'TRBV29*01' , \
        j_gene_name= 'TRBJ2-2*01' , cdr3_aa = 'CASSPTGQLYF')
    # Only the edges are gauranteed to match the real seq shown below, as insertion codons a 
    # unkown and degenerate
    assert r[0:10] == 'tgtgctagcagccccaccgggcagctctacttt'.upper()[0:10]
    assert r[-10:-1] == 'tgtgctagcagccccaccgggcagctctacttt'.upper()[-10:-1] 
    assert r == 'TGTGCTAGCAGTCCTACCGGGCAGCTCTACTTT'

def test_TCRcodon_small_dataframe_beta():
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    df = clone_df_subset[['v_b_gene','j_b_gene', 'cdr3_b_aa','cdr3_b_nucseq']] 
    syn_nucs = df.apply(lambda r: \
        tc.guess_reverse_translation(r['v_b_gene'], r['j_b_gene'], r['cdr3_b_aa'], verbose = False), axis = 1) 
    len_syn = [len(x) for x in syn_nucs]
    len_real =[len(x) for x in df['cdr3_b_nucseq']]
    assert np.all(len_syn == len_real)

def test_TCRcodon_small_dataframe_alpha():
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    df = clone_df_subset[['v_a_gene','j_a_gene', 'cdr3_a_aa','cdr3_a_nucseq']] 
    syn_nucs = df.apply(lambda r: \
        tc.guess_reverse_translation(r['v_a_gene'], r['j_a_gene'], r['cdr3_a_aa'], verbose = False), axis = 1) 
    # Check that synthestic and real seqs are same length
    len_syn = [len(x) for x in syn_nucs]
    len_real =[len(x) for x in df['cdr3_a_nucseq']]
    assert np.all(len_syn == len_real)

def test_TCRcodon_smal_dataframe_alpha_beta_lots():
    """Bigger Example """
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    df = pd.read_csv("tcrdist/test_files_compact/dash.csv") 
    syn_nucs = df.apply(lambda r: \
        tc.guess_reverse_translation( \
        r['v_b_gene'], r['j_b_gene'], r['cdr3_b_aa'],\
        verbose = False), axis = 1) 

    len_syn = [len(x) for x in syn_nucs]
    len_real =[len(x) for x in df['cdr3_b_nucseq']]
    assert np.all(len_syn == len_real)

    syn_nucs = df.apply(lambda r: \
        tc.guess_reverse_translation(\
            r['v_a_gene'], r['j_a_gene'], r['cdr3_a_aa'],\
            verbose = False), axis = 1) 

    # Check that synthestic and real seqs are same length
    len_syn = [len(x) for x in syn_nucs]
    len_real =[len(x) for x in df['cdr3_a_nucseq']]
    assert np.all(len_syn == len_real)

def test_TCRcodon_smal_dataframe_delta_lots():
    tc = TCRcodon(organism = "human", db_file = "gammadelta_db.tsv")
    df = pd.read_csv("tcrdist/test_files_compact/sant.csv")
    # Sant data Doesn't provide J gene so we are handicapped in that regard, for testing we just guess on
    df['j_g_gene'] = 'TRGJ1*01'
    syn_nucs = df.apply(lambda r: \
                        tc.guess_reverse_translation(\
                        r['v_g_gene'], r['j_g_gene'], r['cdr3_g_aa'],\
                        verbose = False), axis = 1) 
   

    # Check that synthestic and real seqs are same length
    len_syn = [len(x) for x in syn_nucs]
    len_real =[3*len(x) for x in df['cdr3_g_aa']]
    assert np.all(len_syn == len_real)

def test_TCRcodon_smal_dataframe_gama_lots():
    tc = TCRcodon(organism = "human", db_file = "gammadelta_db.tsv")
    df = pd.read_csv("tcrdist/test_files_compact/sant.csv")
    df['j_d_gene'] = [tc.get_best_j_gene(aa_seq = x, verbose = False) for x in df['cdr3_d_aa']]
    df = df[df['v_d_gene'].notna()].copy()
    syn_nucs = df.apply(lambda r: \
                        tc.guess_reverse_translation(\
                        r['v_d_gene'], r['j_d_gene'], r['cdr3_d_aa'],\
                        verbose = False), axis = 1) 

    # Check that synthestic and real seqs are same length
    len_syn = [len(x) for x in syn_nucs]
    len_real =[3*len(x) for x in df['cdr3_d_aa']]
    assert np.all(len_syn == len_real)
    assert np.all(len_syn == len_real)

def test_get_bets_j_gene():
    from tcrdist.pairwise import hm_metric                                                                                                                                                                                                                                                    
    tc = TCRcodon(organism = "human", db_file = "gammadelta_db.tsv")
    df = pd.read_csv("tcrdist/test_files_compact/sant.csv")
    someseq = df['cdr3_d_aa'][1]
    x = tc.get_best_j_gene(aa_seq = someseq, verbose = True)
    assert x == 'TRDJ1*01'
    # test for all
    xx = [tc.get_best_j_gene(aa_seq = x, verbose = False) for x in df['cdr3_d_aa']]
    vc = pd.Series(xx).value_counts().to_dict()
    assert vc == {'TRDJ1*01': 271, 'TRDJ3*01': 66, 'TRDJ2*01': 20, 'TRDJ4*01': 9}

def test_get_bets_j_gene_room_for_improvement():
    tc = TCRcodon(organism = "mouse", db_file = "alphabeta_db.tsv")
    df = clone_df_subset[['v_b_gene','j_b_gene', 'cdr3_b_aa','cdr3_b_nucseq']] 
    xx = [tc.get_best_j_gene(aa_seq = x, verbose = False) for x in df['cdr3_b_aa']]
    # THIS CAN'T RESOLVE TIES, OR WHO KNOW WHAT ELSE
    assert np.all(df['j_b_gene'][0:2] == xx[0:2])
    assert np.sum(df['j_b_gene'] == xx) > 20 