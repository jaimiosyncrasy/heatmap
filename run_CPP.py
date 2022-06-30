#region Import statements
import sys
import numpy as np
import my_heatmapSetup_funcs as hm
import my_process_funcs as prc  # new
#endregion

def exp4a(v):
    # ----- good vs. bad branch scenario, 1-step heatmap
    set_acts = ['bus_8','bus_53','bus_57','bus_66'] # has bad branch
    addon_acts= ['bus_152'] # arbitrary

    parmObj=hm.configParms()
    parmObj.set_ctrlTypes(['PBC']*len(set_acts+addon_acts))
    parmObj.set_version(1) # 1 for PBC

    # here we set act_locs as the existing actuators
    feas_configs, heatMapNames=prc.find_good_colocated(parmObj,v.feeder, set_acts, addon_acts, v.node_index_map, v.substation_name, v.depths,v.file_name, v.Vbase_ll, v.Sbase)

    return

def exp4b(v):

    # ----- good vs. bad branch scenario, 1-step heatmap
    set_acts = ['bus_8','bus_53','bus_57','bus_74'] # has good branch
    addon_acts= ['bus_152'] # arbitrary

    parmObj=hm.configParms()
    parmObj.set_ctrlTypes(['PBC']*len(set_acts+addon_acts))
    parmObj.set_version(1) # 1 for PBC

    # here we set act_locs as the existing actuators
    feas_configs, heatMapNames=prc.find_good_colocated(parmObj,v.feeder, set_acts, addon_acts, v.node_index_map, v.substation_name, v.depths,v.file_name, v.Vbase_ll, v.Sbase)

    # result looks like: 'heatmap_colocated' + '_' + file_name

    return