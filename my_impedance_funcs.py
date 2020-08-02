import importlib
import setup_nx # your own module, setup.nx.py
import numpy as np
import math as m
import statistics as st
import cmath
import matplotlib.pyplot as plt 
import itertools
from operator import add
importlib.reload(setup_nx)
from setup_nx import *
from graphviz import Source, render

import datetime
import time

import matplotlib.pyplot as plt
import my_feeder_funcs as ff

#Notes on the following two methods:
#1) Does not account for the presence of artifically added transformer windings, because I would have to separate out
#the line impedance information
#2) Will not function properly if there are cycles
#3) Assumes that switches and transformers have zero impedance; the latter case is due to unit discrepancy: transformer
#impedances use PU, line impedances use Ohms

def get_total_impedance_from_substation(feeder, node_name, depths):
    node_name_full = node_name

    #Notice that since this is a tree, any node will only have at most one predecessor
    total_impedance = np.zeros((3,3), dtype = complex)
    
    current_node = node_name_full
    pred_list = None
    
    try:
        pred_list = list(feeder.network.predecessors(node_name_full))
    except nx.NetworkXError:
        print("Bus with name " + current_node + " does not exist in the feeder.")
        return 0
    
    iter_depth = depths[current_node]
    
    for _ in range(iter_depth):
        #pred_list[0] is the parent node
        impedance = feeder.network.get_edge_data(pred_list[0], current_node, default=None)['connector']
        #print("Type is " + str(type(impedance)))
        if impedance == None:
            print("WARNING: No connection between nodes " + str(pred_list[0]) + " and " + str(current_node) + ".")
            return 0
        else:
            imp_dict = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3), dtype = complex)
            
            total_impedance += imp_dict
            
            current_node = pred_list[0]
            pred_list = list(feeder.network.predecessors(current_node))
    
    return total_impedance

#Method will return a the distance sum of impedances to a common bus upstream if the two buses are not along the
#the same path. For example:
#     A      Calculating total impedance between B and C yields Z_AB + Z_CA
#    / \
#   B   C

def get_total_impedance_between_two_buses(feeder, node_name_1, node_name_2, depths):
    bus_1 = node_name_1
    bus_2 = node_name_2
   
    total_impedance = np.zeros((3,3), dtype = complex)
    
    depth_1 = 0
    depth_2 = 0
    
    try:
        depth_1 = depths[bus_1]
        depth_2 = depths[bus_2]
    except KeyError:
        print("Either the first bus, " + bus_1 + ", or the second bus, " + bus_2 + " is not a valid bus in the feeder.")
        return 0
    
    depth_dif = abs(depth_1 - depth_2)
    
    max_depth_bus = bus_1 if depth_1 > depth_2 else bus_2
    min_depth_bus = bus_1 if max_depth_bus == bus_2 else bus_2
    
    pred_list_max = list(feeder.network.predecessors(max_depth_bus))
    pred_list_min = list(feeder.network.predecessors(min_depth_bus))
    
    for i in range(depth_dif):
        
        impedance = feeder.network.get_edge_data(pred_list_max[0], max_depth_bus, default=None)['connector']
        if impedance == None:
            print("WARNING: No connection between nodes " + str(pred_list_max[0]) + " and " + str(max_depth_bus) + ".")
            return 0
        else:
            imp_dict = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3), dtype = complex)
            
            total_impedance += imp_dict
            
            #Case of where we the two buses are directly linked by purely upstream connections, allowing us to
            #terminate our calculations earlier
            if pred_list_max[0] == min_depth_bus:
                #print("Iterated " + str(i+1) + " times to get direct upstream connection total impedance.")
                return total_impedance
            
            max_depth_bus = pred_list_max[0]
            pred_list_max = list(feeder.network.predecessors(max_depth_bus))
            
    assert(depths[max_depth_bus] == depths[min_depth_bus])
    
    #print("Iterated " + str(depth_dif) + " times to reach equal depths.")
    
    common_parent = pred_list_max[0] == pred_list_min[0]
    
    count_get_to_common = 0
    
    #Here, we simultaneously shift both buses (after the max depth bus has been shifted to be of equal depth to the min
    #depth bus) to a point where the parent bus is shared
    while not common_parent:
        count_get_to_common += 1
        
        impedance_bus_min = feeder.network.get_edge_data(pred_list_min[0], min_depth_bus, default=None)['connector']
        if impedance_bus_min == None:
            print("WARNING: No connection between nodes " + str(pred_list_min[0]) + " and " + str(min_depth_bus) + ".")
            return 0
        else:
            imp_dict_min = impedance_bus_min.Z if isinstance(impedance_bus_min, setup_nx.line) else np.zeros((3,3), dtype = complex)
            
            total_impedance += imp_dict_min
        
            min_depth_bus = pred_list_min[0]
            pred_list_min = list(feeder.network.predecessors(min_depth_bus))
            
        impedance_bus_max = feeder.network.get_edge_data(pred_list_max[0], max_depth_bus, default=None)['connector']
        if impedance_bus_max == None:
            print("WARNING: No connection between nodes " + str(pred_list_max[0]) + " and " + str(max_depth_bus) + ".")
            return 0
        else:
            imp_dict_max = impedance_bus_max.Z if isinstance(impedance_bus_max, setup_nx.line) else np.zeros((3,3), dtype = complex)
            
            total_impedance += imp_dict_max
        
            max_depth_bus = pred_list_max[0]
            pred_list_max = list(feeder.network.predecessors(max_depth_bus))
            
        common_parent = pred_list_max[0] == pred_list_min[0]
    
    print("Total iterations to get to common parent is " + str(count_get_to_common))
    print("Common parent is " + str(pred_list_max[0]))
    
    #Need to iterate one more time to account for "joining" node
    impedance_bus_min = feeder.network.get_edge_data(pred_list_min[0], min_depth_bus, default=None)['connector']
    if impedance_bus_min == None:
        print("WARNING: No connection between nodes " + str(pred_list_min[0]) + " and " + str(min_depth_bus) + ".")
        return 0
    else:
        imp_dict_min = impedance_bus_min.Z if isinstance(impedance_bus_min, setup_nx.line) else np.zeros((3,3), dtype = complex)

        total_impedance += imp_dict_min
    
    impedance_bus_max = feeder.network.get_edge_data(pred_list_max[0], max_depth_bus, default=None)['connector']
    if impedance_bus_max == None:
        print("WARNING: No connection between nodes " + str(pred_list_max[0]) + " and " + str(max_depth_bus) + ".")
        return 0
    else:
        imp_dict_max = impedance_bus_max.Z if isinstance(impedance_bus_max, setup_nx.line) else np.zeros((3,3), dtype = complex)

        total_impedance += imp_dict_max
        
    return total_impedance
    
#Returns the R/X ratio from a node up to the substation       
def get_RX_ratio_tosubst(feeder, node_name,depths):
    impedances_per_phase = get_total_impedance_from_substation(feeder, node_name,depths)
    
    #p1_z = impedances_per_phase['Phase 1']
    #p2_z = impedances_per_phase['Phase 2']
    #p3_z = impedances_per_phase['Phase 3']
    p1_z = impedances_per_phase[0][0]
    p2_z = impedances_per_phase[1][1]
    p3_z = impedances_per_phase[2][2]
    
    p1_rx = np.real(p1_z) / np.imag(p1_z) # R/X ratio
    p2_rx = np.real(p2_z) / np.imag(p1_z) 
    p3_rx = np.real(p3_z) / np.imag(p1_z) 
    
    return {'Phase 1' : np.around(p1_rx,2), 'Phase 2' : np.around(p2_rx,2), 'Phase 3' : np.around(p3_rx,2)}

#Returns the R/X ratio between 2 nodes      
def get_RX_ratio_between_two_buses(feeder, node_name_1, node_name_2,depths):
    impedances_per_phase = get_total_impedance_between_two_buses(feeder, node_name_1, node_name_2,depths)
    
    #p1_z = impedances_per_phase['Phase 1']
    #p2_z = impedances_per_phase['Phase 2']
    #p3_z = impedances_per_phase['Phase 3']
    p1_z = impedances_per_phase[0][0]
    p2_z = impedances_per_phase[1][1]
    p3_z = impedances_per_phase[2][2]

    p1_rx = np.real(p1_z) / np.imag(p1_z) # R/X ratio
    p2_rx = np.real(p2_z) / np.imag(p1_z) 
    p3_rx = np.real(p3_z) / np.imag(p1_z) 
    
    return {'Phase 1' : np.around(p1_rx,2), 'Phase 2' : np.around(p2_rx,2), 'Phase 3' : np.around(p3_rx,2)}

def plot_histogram_RX_ratios(feeder, leaves,slack_bus,depths,leaves_only=False):
    nodes_to_use = leaves if leaves_only else list(feeder.network.nodes)
    
    #PL0001 Case
    #lst_remove = ['bus_N_L_22666_sec', 'bus_N_L_21316_sec', 'bus_N_L_22077_sec', 
                  #'bus_N_L_52586_sec', 'bus_N_L_38426_sec', slack_bus]
    #AL0001 Case
    #lst_remove = ['bus_N_L_87632_sec', 'bus_N_L_46793_sec', 'bus_N_L_17532_sec',
     #            'bus_N_L_108238_sec', 'bus_N_L_111610_sec', 'bus_N_L_147411_sec',
     #            'bus_N_L_6860_sec', 'bus_N_L_9709_sec', 'bus_N_L_110483_sec',
     #            'bus_N_L_132901_sec', 'bus_N_L_40688_sec', 'bus_N_L_138443_sec',
     #            'bus_N_L_108120_sec', 'bus_N_L_116563_sec', 'bus_N_L_113827_sec', slack_bus]
    
    #13NF case
    lst_remove=[]
    
    nodes_to_use = [nd for nd in nodes_to_use if nd not in lst_remove]
    
    RX_ratios_all = [get_RX_ratio_tosubst(feeder, node_name,depths) for node_name in nodes_to_use]
    
    phase_1_RX_ratios = [phase_dct['Phase 1'] for phase_dct in RX_ratios_all]
    phase_2_RX_ratios = [phase_dct['Phase 2'] for phase_dct in RX_ratios_all]
    phase_3_RX_ratios = [phase_dct['Phase 3'] for phase_dct in RX_ratios_all]
    
    plt.title("Phase 1 R/X Ratio Histogram")
    plt.xlabel("R/X Ratio")
    plt.ylabel("Number of Nodes")
    plt.hist(phase_1_RX_ratios, bins=20)
    plt.savefig('AL0001_RXhist_phA.png') # modify for each feeder
    plt.show() # need to savefig before plt.show

    plt.title("Phase 2 R/X Ratio Histogram")
    plt.xlabel("R/X Ratio")
    plt.ylabel("Number of Nodes")
    plt.hist(phase_2_RX_ratios, bins=20)
    plt.savefig('AL0001_RXhist_phB.png') # modify for each feeder
    plt.show()
    
    plt.title("Phase 3 R/X Ratio Histogram")
    plt.xlabel("R/X Ratio")
    plt.ylabel("Number of Nodes")
    plt.hist(phase_3_RX_ratios, bins=20)
    plt.savefig('AL0001_RXhist_phC.png') # modify for each feeder
    plt.show()
    
    return None
