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

import my_feeder_funcs as ff
import my_impedance_funcs as imp
import my_configVis_funcs as vis
import my_detControlMatExistence_funcs as ctrl
import my_detLznRange_funcs as lzn


def get_path_to_substation(feeder, node,depths):
    #returns list of edges, not impedances, between node and substation
    #node = node name as string
    #feeder = initiaized feeder object
    graph = feeder.network
    node_path = []
    current_node = node
    
    for i in range(depths[node]):
        pred_node = list(graph.predecessors(current_node)) #retrieves parent node of current_node
        node_path += [(current_node, pred_node[0])]
        current_node = pred_node[0]
    return node_path


def createRXmatrices(feeder, node_index_map,depths):
    #calculated using phase A impedances unless there is no phase A then phase B else phase C
    #feeder = initiaized feeder object
    #node_index_map = dictionary of node indices with node names as keys
    graph = feeder.network
    n = len(graph.nodes) #number of nodes
    R = np.zeros((n,n)) #initializing R matrix
    X = np.zeros((n,n)) #initializing X matrix
    P = {} #initializing line path dictionary
   
    for node in graph.nodes:
        P[node] = get_path_to_substation(feeder, node,depths)
    
    for n_outer in graph.nodes: #outer loop
        index_outer = node_index_map[n_outer]
        
        for n_inner in graph.nodes: #inner loop
            index_inner = node_index_map[n_inner]
            intersection_set = set.intersection(set(P[n_outer]), set(P[n_inner])) #finds common edges in paths from node to substation
            intersection_list = list(intersection_set)
            r_sum = 0
            x_sum = 0
            
            for edge in intersection_list: #iterates through shared edges and sums their impedances
                tot_edge_impedance = imp.get_total_impedance_between_two_buses(feeder, edge[0], edge[1],depths)
                
                if tot_edge_impedance['Phase 1']:
                    phase_edge_impedance = tot_edge_impedance['Phase 1']
                elif tot_edge_impedance['Phase 2']:
                    phase_edge_impedance = tot_edge_impedance['Phase 2']
                else:
                    phase_edge_impedance = tot_edge_impedance['Phase 3']  
                    
                r_sum += phase_edge_impedance.real
                x_sum += phase_edge_impedance.imag
            R[index_outer][index_inner] = 2*r_sum
            X[index_outer][index_inner] = 2*x_sum
            
    return R, X


def createRXmatrices_3ph(feeder, node_index_map,depths):
    #includes all 3 phases
    #feeder = initiaized feeder object
    #node_index_map = dictionary of node indices with node names as keys
    graph = feeder.network
    n = len(graph.nodes) #number of nodes
    R = np.zeros((3*n, 3*n)) #initializing R matrix
    X = np.zeros((3*n, 3*n)) #initializing X matrix
    P = {} #initializing line path dictionary
   
    for node in graph.nodes:
        P[node] = get_path_to_substation(feeder, node,depths)
    
    for n_outer in graph.nodes: #outer loop
        index_outer = node_index_map[n_outer]
        
        for n_inner in graph.nodes: #inner loop
            index_inner = node_index_map[n_inner]
            intersection_set = set.intersection(set(P[n_outer]), set(P[n_inner])) #finds common edges in paths from node to substation
            intersection_list = list(intersection_set)
            imped_sum = np.zeros((3,3))
            imped_sum = imped_sum.astype('complex128')
            
            for edge in intersection_list: #iterates through shared edges and sums their impedances
                impedance = feeder.network.get_edge_data(edge[1], edge[0], default=None)['connector']
                tot_edge_impedance = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3)) 
                imped_sum += tot_edge_impedance
            #print('impedance.Z=',impedance.Z)
            for i_row in range(3):    
                for i_col in range(3):
                    R[(3 * index_outer) + i_row][(3 * index_inner) + i_col] = 2 * imped_sum[i_row][i_col].real
                    X[(3 * index_outer) + i_row][(3 * index_inner) + i_col] = 2 * imped_sum[i_row][i_col].imag
            
    return R, X


def createNodeIndexMap(feeder):
    #indexes assigned based on numerical value of node, not node's possition in the network
    #for example, the node named with the lowest numerical value will map to 0, not the substation
    #feeder = initiaized feeder object
    graph = feeder.network
    node_index_map = {} #node indices for indicMat and F matrix
    t = 0 #initializing first index for node_index_map
    
    for node in graph.nodes: #populating node_idex_map
        node_index_map[node] = t
        t += 1
    return node_index_map


def getKey(dictionary, value):
    #use to retrieve node name from node_index_map using node index
    #value = node index
    for k, v in dictionary.items(): 
         if value == v:
                key = k
    return key 


def setupStateSpace(n, feeder, node_index_map,depths):
    #n = number of nodes in network
    #feeder = initiaized feeder object
    #node_index_map = dictionary of node indices with node names as keys
    A = np.identity(6*n)
    R, X = createRXmatrices_3ph(feeder, node_index_map,depths)
    concat_XR = np.concatenate((X, R), axis = 1)
    concat_XR_halfs = np.concatenate(((-1/2) * R, (1/2) * X), axis = 1)
    B = np.concatenate((concat_XR, concat_XR_halfs))
    return A, B


# correct version
def computeFeas_v1(feeder, act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map):
    node_1 = list(feeder.network.successors(substation_name))
    z12 = imp.get_total_impedance_from_substation(feeder, node_1[0],depths) # 3 phase
    MYfeas,MYfeasFs,MYnumfeas,MYnumTried,MYnumact = ctrl.detControlMatExistence(feeder, act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map)
    print('num feas=',MYnumfeas)
    print('num tried=',MYnumTried)

    #maxLznError = lzn.detLznRange(feeder, Vbase_ll, Sbase, z12, act_locs, load_data, headerpath, substation_name)
    #return matExist, maxLznError
    MYmaxerror=-1 # temporary
    return MYfeas,MYmaxerror

# workaround version
def computeFeas_v2(feeder, act_locs, perf_nodes, A, B, indicMat):
    random_bit = np.random.choice(2,1)
    if random_bit[0] == 0:
        feas = True
    else:
        feas = False
    maxError = 10
    return feas, maxError


def updateStateSpace(n, act_locs, perf_nodes, node_index_map):
    #creates (2n*2n) matrix with 1 and ij and (i+n)(j+n) if there is an actuator at i tracking a performance node at j, and 0 otherwise
    #n = number of nodes in network
    #act_locs = list of actuators locations in network (list of strings)
    #perf_node = performance node as string
    #node_index_map = dictionary of node indices for indicMat and F matrix
    indicMat = np.zeros((6*n,6*n))
    
    for i in range(len(act_locs)): 
        act = act_locs[i]
        perf = perf_nodes[i]
        act_index = node_index_map[act]
        perf_index = node_index_map[perf]
        
 #       for i_row in range(3):    
 #           for i_col in range(3):
 #               indicMat[act_index*3+i_row][perf_index*3+i_col] = 1
 #               indicMat[act_index*3+3*n+i_row][perf_index*3+3*n+i_col] = 1                  
        for k in range(3):
            indicMat[act_index*3+k][perf_index*3+k] = 1
            indicMat[act_index*3+3*n+k][perf_index*3+3*n+k] = 1   
            
    return indicMat
    
#     for i in range(len(act_locs)): 
#         act = act_locs[i]
#         perf = perf_nodes[i]
#         act_index = node_index_map[act]
#         perf_index = node_index_map[perf]
#         indicMat[act_index][perf_index] = 1
#         indicMat[act_index+n][perf_index+n] = 1
#     return indicMat
    
def markActLoc(graph, act_loc):
    #changes color of nodes with set actuators to turquoise
    #graph = networkx graph object (feeder.network)
    #act_loc = node name as string where actuator is placed
    graph.nodes[act_loc]['style'] = 'filled'
    graph.nodes[act_loc]['fillcolor'] = 'turquoise'
    return
        

def markFeas(feas, test_act_loc, graph):
    #if controllability can be achieved with actuator at test_act_loc then mark green, otherwise mark red
    #feas = True or False
    #test_act_loc = node name as string
    #graph = networkx graph object (feeder.network)
    graph.nodes[test_act_loc]['style'] = 'filled'
   
    if feas:
        graph.nodes[test_act_loc]['fillcolor'] = 'green'
    else:
        graph.nodes[test_act_loc]['fillcolor'] = 'red'
    return

def eval_config(feeder, all_act_locs, perf_nodes, node_index_map,substation_name,depths,file_name):
    #all_act_locs and perf_nodes = lists of node names as strings
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map,depths)
    indicMat = updateStateSpace(n, all_act_locs, perf_nodes, node_index_map)
    feas, maxError = computeFeas_v1(feeder, all_act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map)
    vis.markActuatorConfig(all_act_locs, feeder, file_name) # create diagram with actuator locs marked
    
    print('Actuator configuration is feasible') if feas else print('Actuator configuration is not feasible')
    return feas, maxError

def find_good_colocated(feeder, node_index_map,substation_name,depths, file_name):
    # almost the same as runheatmap process, but only runs once and shows one heatmap indicating which nodes are good to place a co-located act/perf node
    # return list of all "green" configs and the associated lzn errors
    
    #all_act_locs and perf_nodes = lists of node names as strings
    a = 0
    graph = feeder.network
    heatMapNames=[] # collect heat map names as list of strings
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map,depths)
    lzn_error_run_sum = 0
    feas_configs=[] 
    lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
    test_nodes = []
    graphNodes_nosub=np.delete(graph.nodes,[6, 7]) # dont consider co-located at substation nodes, node 650 and 651
    # Note: ^idx 6 & 7 are MANUALLY PICKED OUT FOR 13NF
    
    for node in graphNodes_nosub: # try placing act/perf at all nodes of the network
        test_nodes.append(node)

    for test in test_nodes:
        indicMat = updateStateSpace(n, [test], [test], node_index_map) # (n,act,perf,dictionary)
        print('evaluating act at ',[test],', perf at ',[test])
        feas, maxError = computeFeas_v1(feeder, [test], A, B, indicMat,substation_name,[test],depths,node_index_map) # pass in potential actual loc
        lzn_error_dic[test] = maxError
        markFeas(feas, test, graph)
        if feas:
            feas_dic = {}
            feas_dic['act'] = [test]
            feas_dic['perf'] = [test]
            feas_dic['lznErr'] = [lzn_error_dic[test]]
            feas_configs += [feas_dic]        

    heatMapName='heat_map_colocated' + '_' + file_name
    heatMapNames.append(heatMapName)
    nx.nx_pydot.write_dot(graph, heatMapName)
    render('dot', 'png', heatMapName)

    return feas_configs, heatMapNames

def runHeatMapProcess(feeder, all_act_locs, perf_nodes, node_index_map,substation_name,depths, file_name):
    # On empty network, compute heatmap (assess feas and lzn error on every node of the feeder, color each red/green on diagram)
    # Then for each act-perf node pair, compute heatmap
    # return list of all "green" configs and the associated lzn errors
    # heatmap shows good places to placed act wrt given perf node, NOT good places to place colocated act-perf node
    
    #all_act_locs and perf_nodes = lists of node names as strings
    a = 0
    graph = feeder.network
    cur_act_locs = []
    cur_perf_nodes=[]
    heatMapNames=[] # collect heat map names as list of strings
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map,depths)
    lzn_error_run_sum = 0
    feas_configs=[]
    
    while a < len(all_act_locs): #outer loop, a=number of actuators
        for act in cur_act_locs: 
            markActLoc(graph, act)
            
        lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
        test_nodes = []
        graphNodes_nosub=np.delete(graph.nodes,[6, 7]) # dont consider co-located at substation nodes, node 650 and 651 
        for node in graphNodes_nosub:
            if node not in cur_act_locs:
                test_nodes.append(node)
    
        for test in test_nodes: #inner loop
            # heatmap color indicates indicates good places to place actuator given chosen loc of perf node (not necessarily colocated)
            indicMat = updateStateSpace(n, [test] + cur_act_locs, [perf_nodes[a]] + cur_perf_nodes, node_index_map)
            print('evaluating act at ',[test]+cur_act_locs,', perf at ',[perf_nodes[a]] + cur_perf_nodes)

            feas, maxError = computeFeas_v1(feeder, [test]+cur_act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map)
            lzn_error_dic[test] = maxError
            markFeas(feas, test, graph)
            if feas:
                feas_dic = {}
                feas_dic['act'] = [test] + cur_act_locs
                feas_dic['perf'] = [test] + cur_perf_nodes
                feas_dic['lznErr'] = [lzn_error_dic[test]]
                feas_configs += [feas_dic]        
        
        heatMapName='heat_map_' + str(a) + '_' + file_name
        heatMapNames.append(heatMapName)
        nx.nx_pydot.write_dot(graph, heatMapName)
        render('dot', 'png', heatMapName)
        a += 1
        
        if a <= len(all_act_locs): # choose actuator and finalize assoc perf node
            cur_act_locs = all_act_locs[0:a] # populate cur_act_locs with subset of all_act_locs
            cur_perf_nodes = perf_nodes[0:a] 
            lzn_error_run_sum += lzn_error_dic[cur_act_locs[-1]]
            print('The total max linearization error after '+ str(a) +' actuators have been placed = '+ str(lzn_error_run_sum))
            
        # end of while loop
    return feas_configs, lzn_error_run_sum, heatMapNames

