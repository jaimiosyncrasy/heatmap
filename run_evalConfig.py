#region Import statements
import sys
import numpy as np
import my_impedance_funcs as imp
import my_heatmapSetup_funcs as hm
import my_process_funcs as prc  # new
import my_configVis_funcs as vis
import pickle
import my_feeder_funcs as ff

#endregion

#exp 1
cfg1 = ['bus_300', 'bus_197', 'bus_104', 'bus_108', 'bus_106']  # bad config 1, cluster in one portion of feeder
cfg2 = ['bus_76', 'bus_49', 'bus_109', 'bus_251',
        'bus_66']  # good config 1, same z2sub as config 1 but greater stability
cfg3 = ['bus_10', 'bus_85', 'bus_56', 'bus_350', 'bus_43']  # bad config 2, evenly spaced

# exp2
cfg5 = ['bus_87', 'bus_90', 'bus_92', 'bus_95', 'bus_80', 'bus_73', 'bus_70', 'bus_67', 'bus_160', 'bus_61']


def exp1a(v): # run bad config 1 and good config 1

    # --------- section 1.a: run eval_config on CPP config --------------
    parmObj = hm.configParms()
    parmObj.set_version(1)  # 1 for PBC

    # compute sum of z2subs for a given config, show they are about the same
    z2sub_cfg1 = imp.get_z2sub_for_config(v.feeder, cfg1, v.depths)
    print(np.around(z2sub_cfg1, 2))
    z2sub_cfg2 = imp.get_z2sub_for_config(v.feeder, cfg2, v.depths)
    print(np.around(z2sub_cfg2, 2))

    parmObj.set_ctrlTypes(['PBC'] * len(cfg1))
    feas, maxError, percent_feas, bestF_asvec,bestF_asmat, indicMat = prc.eval_config(parmObj, v.feeder, cfg1, cfg1, v.node_index_map,
                                                               v.substation_name, v.depths, v.file_name, v.Vbase_ll, v.Sbase)
    # Saving results:
    with open('exp1_cfg1_7.3.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([feas, percent_feas, bestF_asvec,bestF_asmat], f)


    parmObj.set_ctrlTypes(['PBC'] * len(cfg2))
    feas, maxError, percent_feas, bestF_asvec,bestF_asmat, indicMat = prc.eval_config(parmObj, v.feeder, cfg2, cfg2, v.node_index_map,
                                                               v.substation_name, v.depths, v.file_name, v.Vbase_ll, v.Sbase)
    # Saving results:
    with open('exp1_cfg2_7.3.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([feas, percent_feas, bestF_asvec,bestF_asmat], f)


    return

def exp1b(v): # run bad config 2

    # --------- section 1.b: run eval_config on CPP config --------------
    parmObj = hm.configParms()
    parmObj.set_version(1)  # 1 for PBC

    # cfg1=['bus_18','bus_26','bus_38','bus_58','bus_63','bus_73','bus_99','bus_91','bus_109'] # bad config 1
    # cfg2=['bus_18','bus_26','bus_38','bus_58','bus_63','bus_73','bus_99','bus_91','bus_109'] # good config 1

    # # compute sum of z2subs for a given config
    # z2sub_cfg1=imp.get_z2sub_for_config(feeder,cfg1,depths)
    # print(np.around(z2sub_cfg1,2))
    # z2sub_cfg2=imp.get_z2sub_for_config(feeder,cfg2,depths)
    # print(np.around(z2sub_cfg2,2))

    parmObj.set_ctrlTypes(['PBC'] * len(cfg3))
    feas, maxError, percent_feas, bestF_asvec,bestF_asmat, indicMat = prc.eval_config(parmObj, v.feeder, cfg3, cfg3, v.node_index_map,
                                                               v.substation_name, v.depths, v.file_name, v.Vbase_ll, v.Sbase)
    # Saving results:
    with open('exp1_cfg3_7.3.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([feas, percent_feas, bestF_asvec,bestF_asmat], f)
    return

def exp2(v):
    # --------- section 3: run eval_config on 1-sensor config --------------
    # Need to design non-colocated configs:

    # 1 sensor deep in feeder, 1 DER on top of it
    # 10 DERs nearby, some above and some below
    # DERs all non-sensor nodes every other node

    parmObj = hm.configParms()
    parmObj.set_version(1)  # 1 for PBC

    perf_node = ['bus_76']
    cfg4 = perf_node
    cfg6 = v.graphNodes_nosub
    cfg6.remove(perf_node[0])

    parmObj.set_ctrlTypes(['PBC'] * len(cfg4))
    feas, maxError, percent_feas, bestF_asvec,bestF_asmat, indicMat = prc.eval_config(parmObj, v.feeder, cfg4, perf_node*len(cfg4), v.node_index_map,
                                                               v.substation_name, v.depths, v.file_name, v.Vbase_ll, v.Sbase)
    with open('exp2_cfg4_7.3.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([feas, percent_feas, bestF_asvec,bestF_asmat], f)
    parmObj.set_ctrlTypes(['PBC'] * len(cfg5))
    feas, maxError, percent_feas, bestF_asvec,bestF_asmat, indicMat = prc.eval_config(parmObj, v.feeder, cfg5, perf_node*len(cfg5), v.node_index_map,
                                                               v.substation_name, v.depths, v.file_name, v.Vbase_ll, v.Sbase)
    with open('exp2_cfg5_7.3.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([feas, percent_feas, bestF_asvec,bestF_asmat], f)
    parmObj.set_ctrlTypes(['PBC'] * len(cfg6))
    feas, maxError, percent_feas, bestF_asvec,bestF_asmat, indicMat = prc.eval_config(parmObj, v.feeder, cfg6, perf_node*len(cfg6), v.node_index_map,
                                                               v.substation_name, v.depths, v.file_name, v.Vbase_ll, v.Sbase)
    with open('exp2_cfg6_7.3.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([feas, percent_feas, bestF_asvec,bestF_asmat], f)
    return

def exp1_markGraph(v):

    with open('exp1_cfg2_7.3.pkl','rb') as f:  #open results from the good config
        [feas, percent_feas, bestF_asvec, bestF_asmat] = pickle.load(f)

    vis.markMultiple_actConfig([cfg1,cfg2,cfg3], v.feeder, v.file_name)
    return

def exp2_markGraph(v):
    ff.clear_graph(v.feeder)
    vis.markActuatorConfig(cfg5, v.feeder, v.file_name)  # create diagram with actuator locs marked
    return