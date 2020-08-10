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
import my_impedance_funcs as imp
import my_detControlMatExistence_funcs as ctrl
import my_detLznRange_funcs as lzn
import my_heatmapSetup_funcs as hm


def markActuatorConfig(lst_act_locs, feeder, file_name):
    # creates diagram of actuator-configuration (though leaves out performance nodes)
    # lst_act_locs = list of node names where actuators are placed
    # feeder = initialized feeder object
    # file_name = string that will be used as the file name of the returned network graph
    ff.clear_graph(feeder)
    graph = feeder.network
    
    for loc in lst_act_locs:
        graph.nodes[loc]['style'] = 'filled'
        graph.nodes[loc]['fillcolor'] = 'gray'
    nx.nx_pydot.write_dot(graph, 'actConfig_'+ file_name)
    render('dot', 'png', 'actConfig_'+ file_name)
    return


def markCommonFeasNodes(lst_feas_configs, feeder):
    # input a list of feasible actuator configurations where each configuration is its own list (ie a list of lists)
    # for actuator locations that appear in every configuration, mark node blue
    shared_locs = lst_feas_configs[0]
    num_configs = len(lst_feas_configs)
    graph = feeder.network
    
    for i in range(num_configs):
        if i <= num_configs - 2:
            cur_config = lst_feas_configs[i + 1]
            copy_shared_locs = shared_locs
            for loc in copy_shared_locs:
                if loc not in cur_config:
                    shared_locs.remove(loc)
    
    for act_loc in shared_locs:
        graph.nodes[act_loc]['style'] = 'filled'
        graph.nodes[act_loc]['fillcolor'] = 'gray'
        
    nx.nx_pydot.write_dot(graph, 'shared_act_locs_' + file_name)
    render('dot', 'png', 'shared_act_locs_' + file_name)
    return shared_locs


def find_ActLoc_in_FeasConfigs(act_loc, feas_configs):
    # act_loc = node name of potential actuator location
    # feas_configs = a list of feasible actuator configurations
    configs_with_actLoc = []
    
    for config in feas_configs:
        if act_loc in config:
            configs_with_actLoc += [config]
    return configs_with_actLoc
    

def plot_actuator_num_histogram(lst_feas_configs):
    # input list of lists of feasible actuator configurations
    # output histogram reflecting number of actuators in each configuration
    lst_num_acts = []
    for config in lst_feas_configs:
        lst_num_acts += [len(config)]
        
    bins_lst = range(max(lst_num_acts) + 2)
        
    plt.title("Number of Actuators per Configuration")
    plt.xlabel("Number of Actuators per Configuration")
    plt.ylabel("Number of Feasible Configurations")
    plt.hist(lst_num_acts, bins = bins_lst, edgecolor = 'black', align = 'left')
    plt.savefig(file_name + '_' + 'Histogram.png') 
    plt.show() # need to savefig before plt.show
    return


#used in multiphase line loss func
def assign_network_branches1(feeder, substation_name):
    # feeder = initialized feeder
    # substation_name = name of substation node as string
    cur_child_nodes = list(feeder.network.successors(substation_name))
    branches = []
    branch_builder = [substation_name]
    
    while len(cur_child_nodes) == 1:
        branch_builder += cur_child_nodes
        cur_child_nodes = list(feeder.network.successors(cur_child_nodes[0]))
            
    if len(cur_child_nodes) > 1:
        branches += [branch_builder]
        for child in cur_child_nodes:
            branches += assign_network_branches1(feeder, child)
            
    else:
        branches += [branch_builder]
        return branches
    
    return branches  


# used for general visualization
def assign_network_branches(feeder, substation_name):
    # feeder = initialized feeder
    # substation_name = name of substation node as string
    cur_child_node = substation_name
    branches = []
    branch_builder = []
    
    while cur_child_node:
        branch_builder += [cur_child_node]
        all_children = list(feeder.network.successors(cur_child_node))
        
        if len(all_children) > 1:
            cur_child_node = max(all_children)
            all_children.remove(cur_child_node)
            new_branch_heads = all_children
            
            for head in new_branch_heads:
                branches += assign_network_branches3(feeder, head)
                    
        elif len(all_children) == 1:
            cur_child_node = all_children[0]
            
        else:
            break
    
    branches += [branch_builder]
    return branches


def mark_network_branches(feeder, branch_lst):
    # feeder = initialized feeder object
    # branch_lst = a list of lists containing node names as strings of nodes that are in the same branch
    # can mark up to 21 uniquely colored branches, for branch lists with more than 21 branches colors will repeat
    graph = feeder.network
    cur_color = 1
    color_scheme = 'set312'
    
    for branch in branch_lst:
        
        for node in branch:
            graph.nodes[node]['colorscheme'] = color_scheme
            graph.nodes[node]['style'] = 'filled'
            graph.nodes[node]['fillcolor'] = cur_color
        
        cur_color += 1
        
        if color_scheme == 'set312' and cur_color > 12:
            cur_color = 1
            color_scheme = 'set19'
        elif color_scheme == 'set19' and cur_color > 9:
            cur_color = 1
            color_scheme = 'set312'
    
    nx.nx_pydot.write_dot(graph, 'branch_key:' + file_name)
    render('dot', 'png', 'branch_key:' + file_name)
    return
        

def find_good_branches(lst_feas_configs, branch_lst, num_good_branches):
    # branch_lst = list of branches for a certain feeder --> get by calling assign_network_branches()
    # lst_feas_configs = list of feasible actuator configurations --> get by calling runHeatMapProcess()
    # num_good_branches = number of 'good' branches you want returned by the function, argument should be an integer, branches returned in order from 'best' to 'worst'
    # the 'best' branch is the branch with the most actuators on it across all configurations
    branch_dic_unique = {}
    branch_dic_repeat = {}
    best_branches_unique = []
    best_branches_repeat = []
    
    for branch in branch_lst:
        branch_head = branch[0]
        branch_dic_unique[branch_head] = 0
        branch_dic_repeat[branch_head] = 0
    
    for feas in lst_feas_configs:
        unique_tracker = 0 # tracks whether or not an a branch has already been represented in a certain configuration
       
        for act_loc in feas:
            
            for branch in branch_lst:
                branch_head = branch[0]
                
                if act_loc in branch and unique_tracker == 0:
                    branch_dic_unique[branch_head] += 1
                    branch_dic_repeat[branch_head] += 1
                    unique_tracker = 1
                    break
                    
                elif act_loc in branch and unique_tracker == 1:
                    branch_dic_repeat[branch_head] += 1
    
    common_branches = [k for k, v in branch_dic_unique.items() if v > 0]
    counter = num_good_branches
    while counter > 0:
        max_unique = max(branch_dic_unique, key = branch_dic_unique.get)
        max_repeat = max(branch_dic_repeat, key = branch_dic_repeat.get)
        best_branches_unique += [max_unique + ' = ' + str(branch_dic_unique[max_unique]) + ' actuators']
        best_branches_repeat += [max_repeat + ' = ' + str(branch_dic_repeat[max_repeat]) + ' actuators']
        branch_dic_unique.pop(max_unique)
        branch_dic_repeat.pop(max_repeat)
        counter -= 1
    
    print('Branches are represented by their first node (the node closest to the substation).')
    print('\nCommon Branches:')
    print(common_branches)
    print('\nThe ' + str(num_good_branches) + ' best branches when each actuator configuration is only considered once:')
    print(best_branches_unique)
    print('\nThe ' + str(num_good_branches) + ' best branches when each actuator configuration is considered multiple times:')            
    print(best_branches_repeat)
    return


def determine_good_or_bad_branch(branch, lst_feas_configs, num_configs_for_good_branch):
    # branch = list of nodes in a branch
    # lst_feas_configs = list of feasible actuator configurations
    # num_configs_for_good_branch = the minimum number of configurations that must use a branch for the branch to be 'good' (argument should be an integer)
    # a configuration 'uses' a branch if the configuration includes an actuator that is placed on one of the nodes within the branch
    configs_with_branch = []
    
    for config in lst_feas_configs:
        
        for node in config:
            
            if node in branch:
                configs_with_branch += [config]
                break
                
    if len(configs_with_branch) >= num_configs_for_good_branch:
        print('Branch ' + str(branch) + ' is good.')
        
    else:
        print('Branch ' + str(branch) + ' is bad.')
    
    print('\nNumber of configurations that use the branch: ' + str(len(configs_with_branch)))
    return configs_with_branch          


def find_branch_in_branch_list(node_in_branch, branch_lst):
    # node_in_branch = a node in the branch you want to select
    # branch_lst = list of branches in a network
    for branch in branch_lst:
        
        if node_in_branch in branch:
            my_branch = branch
            break
    return my_branch
            
    
def phaseCouplingPerNode(feeder,depths):
    # phase coupling = mutual impedance/ self impedance
    graph = feeder.network
    nodes = graph.nodes
    coupling_ratios = {}
    
    for node in nodes:
        edge_path = hm.get_path_to_substation(feeder, node,depths)
        self_imped = 0
        mutual_imped = 0
            
        for edge in edge_path:
            impedance_test = graph.get_edge_data(edge[1], edge[0], default=None)['connector']
            impedance = impedance_test.Z if isinstance(impedance_test, setup_nx.line) else np.zeros((3,3))
            self_imped += impedance[0][0] + impedance[1][1] + impedance[2][2]
            mutual_imped += impedance[0][1] + impedance[0][2] + impedance[1][0] + impedance[1][2] + impedance[2][0] + impedance[2][1]
        
        coupling_ratios[node] = mutual_imped / self_imped if self_imped != 0 else mutual_imped
    return coupling_ratios


def createColorMap(feeder, values_dic, file_name):
    graph = feeder.network
    vals = values_dic.values()
    if isinstance(list(vals)[0], complex):
        vals = [(val.imag**2 + val.real**2)**(1/2) for val in vals]
    maxVal = max(vals)
    minVal = min(vals)
    bin_size = (maxVal - minVal)/8
    bin_edges = [[minVal + (i*bin_size), minVal + ((i + 1)*bin_size)] for i in range(8)]
    # color reference: https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib
    colors = ['gray','green', 'limegreen', 'yellowgreen', 'yellow','gold', 'orange', 'red']
    bins_with_clrs = [bin_edges[i] + [colors[i]] for i in range(8)]
    print(bins_with_clrs)
    
    for node, val in values_dic.items():
        if isinstance(val, complex):
            color_num = (val.imag**2 + val.real**2)**(1/2)
        else:
            color_num = val
            
        for b in bins_with_clrs:
            if color_num == minVal and b[0] == minVal:
                graph.nodes[node]['style'] = 'filled'
                graph.nodes[node]['fillcolor'] = b[2]
                break
            elif color_num > b[0] and color_num <= b[1]:
                graph.nodes[node]['style'] = 'filled'
                graph.nodes[node]['fillcolor'] = b[2]
                break
  
    nx.nx_pydot.write_dot(graph, file_name)
    render('dot', 'png', file_name)

    
def createColorMap_BandW(feeder, values_dic, file_name):
    graph = feeder.network
    vals = values_dic.values()
    if isinstance(list(vals)[0], complex):
        vals = [(val.imag**2 + val.real**2)**(1/2) for val in vals]
    maxVal = max(vals)
    minVal = min(vals)
    bin_size = (maxVal - minVal)/8
    bin_edges = [[minVal + (i*bin_size), minVal + ((i + 1)*bin_size)] for i in range(8)]
    sizes = [.1, .2, .3, .4, .5, .6, .7, .8]
    bins_with_sizes = [bin_edges[i] + [sizes[i]] for i in range(8)]
    
    for node, val in values_dic.items():
        if isinstance(val, complex):
            size_num = (val.imag**2 + val.real**2)**(1/2)
        else:
            size_num = val
            
        for b in bins_with_sizes:
            if size_num == minVal and b[0] == minVal:
                graph.nodes[node]['width'] = b[2]
                break
            elif size_num > b[0] and size_num <= b[1]:
                graph.nodes[node]['width'] = b[2]
                break
                
    nx.nx_pydot.write_dot(graph, file_name)