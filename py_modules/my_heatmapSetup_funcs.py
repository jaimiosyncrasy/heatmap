import importlib
import setup_nx # your own module, setup.nx.py
import numpy as np
import math as m
import statistics as st
import cmath
import matplotlib.pyplot as plt 
import itertools
import random
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


def get_path_to_substation(feeder, node, depths):
    #returns list of edges (not impedances) between node and substation
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


def createRXmatrices_3ph(feeder, node_index_map, depths):
    #returns 2 (3n)x(3n) matrices containing R and X values for the network 
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
                impedance = feeder.network.get_edge_data(edge[1], edge[0], default = None)['connector']
                tot_edge_impedance = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3)) 
                imped_sum += tot_edge_impedance
            #print('impedance.Z=',impedance.Z)
            for i_row in range(3):    
                for i_col in range(3):
                    R[(3*index_outer) + i_row][(3*index_inner) + i_col] = 2*imped_sum[i_row][i_col].real
                    X[(3*index_outer) + i_row][(3*index_inner) + i_col] = 2*imped_sum[i_row][i_col].imag
            
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


def setupStateSpace(n, feeder, node_index_map, depths):
    #initializes state space matrices A and B
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
def computeFeas_v1(feeder, act_locs, A, B, indicMat, substation_name, perf_nodes, depths, node_index_map, Vbase_ll, Sbase, load_data, headerpath, modelpath, printCurves):
    node_0 = list(feeder.network.successors(substation_name))
    node_1 = list(feeder.network.successors(node_0[0]))
    z12 = imp.get_total_impedance_from_substation(feeder, node_1[0],depths) # 3 phase, not pu
    B12=np.zeros((3,3)) # TEMPORARY, line susceptance, Yshunt=G+jB

    MYfeas,MYfeasFs,MYnumfeas,MYnumTried,MYnumact = ctrl.detControlMatExistence(feeder, act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map)
    print('num feas=',MYnumfeas)
    print('num tried=',MYnumTried)

    #lzn_err_max, slopes = lzn.detLznRange(feeder, Vbase_ll, Sbase, z12, B12, act_locs, load_data, headerpath, substation_name, modelpath, depths,printCurves) # usually called by computeFeas
    lzn_err_max=[-1, -1, -1, -1] # workaround, for [PV, QV, Pdel,Qdel] lzn errors

    return MYfeas,lzn_err_max,MYnumfeas

# workaround version
def computeFeas_v2(feeder, act_locs, perf_nodes, A, B, indicMat):
    random_bit = np.random.choice(2,1)
    if random_bit[0] == 0:
        feas = True
    else:
        feas = False
    maxError = 10
    return feas, maxError


def updateStateSpace(feeder, n, act_locs, perf_nodes, node_index_map):
    #creates (6n*6n) matrix with 1 at (3i+ph)(3j+ph) and (3i+3n+ph)(3j+3n+ph) if there is an actuator at index i tracking a performance node at index j, and 0 otherwise
    #in the above description, ph is the integer representation (a=0, b=1, c=2) of the phase intersection between the actuator and performance nodes
    #if an actuator and performance node have no phases in common, a warning is printed
    #n = number of nodes in network
    #act_locs = list of actuators locations in network (list of strings)
    #perf_nodes = list of performance nodes 
    #node_index_map = dictionary of node indices for indicMat and F matrix
    indicMat = np.zeros((6*n,6*n))
    for i in range(len(act_locs)): 
        act = act_locs[i]
        perf = perf_nodes[i]
        act_phases = feeder.busdict[act[4:]].phases
        perf_phases = feeder.busdict[perf[4:]].phases
        act_index = node_index_map[act]
        perf_index = node_index_map[perf]
        
        phase_intrsct = [ph for ph in act_phases if ph in perf_phases]
        if phase_intrsct == []:
            print('WARNING: act_node ' + act + ' can NOT track perf_node ' + perf + ' --> no common phases')
        
        for i in range(len(phase_intrsct)):
            if phase_intrsct[i] == 'a':
                phase_intrsct[i] = 0
            elif phase_intrsct[i] == 'b':
                phase_intrsct[i] = 1
            elif phase_intrsct[i] == 'c':
                phase_intrsct[i] = 2
                          
        for ph in phase_intrsct:
            indicMat[(act_index*3)+ph][(perf_index*3)+ph] = 1
            indicMat[(act_index*3)+(3*n)+ph][(perf_index*3)+(3*n)+ph] = 1   
            
    return indicMat

    
def markActLoc(graph, act_loc):
    #changes color of nodes with set actuators to gray
    #graph = networkx graph object (feeder.network)
    #act_loc = node name as string where actuator is placed
    graph.nodes[act_loc]['style'] = 'filled'
    graph.nodes[act_loc]['fillcolor'] = 'gray'
    graph.nodes[act_loc]['shape'] = 'circle'
    return
        

def markFeas(numfeas, test_act_loc, graph):
    #if controllability can be achieved with actuator at test_act_loc then mark green, if only a few controllable configurations (given by thresh_yellowgreen) exist mark yellow, otherwise mark red
    #feas = True or False
    #test_act_loc = node name as string
    #graph = networkx graph object (feeder.network)
    graph.nodes[test_act_loc]['style'] = 'filled'
    graph.nodes[test_act_loc]['shape'] = 'circle'

    thresh_yellowgreen = 15 # you choose
    # for more colors see https://graphviz.org/doc/info/attrs.html
    if numfeas >= thresh_yellowgreen:
        graph.nodes[test_act_loc]['fillcolor'] = 'turquoise'
    elif numfeas >= 1:
        graph.nodes[test_act_loc]['fillcolor'] = 'yellow'
    else:
        graph.nodes[test_act_loc]['fillcolor'] = 'red'
    return


def eval_config(feeder, all_act_locs, perf_nodes, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #returns whether controllability can be achieved for a given actuator performance node configuration
    #also returns the linearization error associated with the feasibility calculation
    #all_act_locs and perf_nodes = lists of node names as strings
    printCurves = True # your choice on whether to print PVcurves
    
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map,depths)
    indicMat = updateStateSpace(feeder, n, all_act_locs, perf_nodes, node_index_map)
    feas, maxError, numfeas = computeFeas_v1(feeder, all_act_locs, A, B, indicMat, substation_name, perf_nodes, depths, node_index_map, Vbase_ll, Sbase, load_data, headerpath, modelpath, printCurves)
    vis.markActuatorConfig(all_act_locs, feeder, file_name) # create diagram with actuator locs marked
    
    print('Actuator configuration is feasible') if feas else print('Actuator configuration is not feasible')
    return feas, maxError, numfeas


def remove_subst_nodes(feeder, file_name):
    #remove the substation nodes/node from the network's node list
    graph = feeder.network
    if file_name == '13NF':
        substIdx = [6, 7] # substation index --> Note: idx 6 & 7 are MANUALLY PICKED OUT FOR 13NF
    elif file_name == '123NF':
        substIdx = [22, 24]
    elif file_name == 'PL0001':
        substIdx = [340] 
    
    graphNodes_nosub = np.delete(graph.nodes, substIdx) # dont consider substation nodes, node 650 and 651 for 13NF
    return graphNodes_nosub

    
def find_good_colocated(feeder, act_locs, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #almost the same as runheatmap process, but only runs once and shows one heatmap indicating which nodes are good to place a co-located act/perf node
    #return list of all "green" configs and the associated lzn errors
    #act_locs == a list of pre-set colocated act/perf locations --> to evaluate an empty network, pass in act_locs == []
    a = 0
    ff.clear_graph(feeder) # clear any previous modifictions made to graph
    graph = feeder.network
    heatMapNames = [] # collect heat map names as list of strings
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map, depths)
    feas_configs = [] 
    lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
    test_nodes = []
    printCurves=False # your choice on whether to print PVcurves
    graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider co-located at substation nodes, node 650 and 651
    
    for node in graphNodes_nosub: # try placing act/perf at all nodes of the network
        if node not in act_locs:
            test_nodes.append(node)
    
    for act in act_locs: 
        markActLoc(graph, act)

    for test in test_nodes:
        print('evaluating act and perf colocated at ',[test]) 
        indicMat = updateStateSpace(feeder, n, [test] + act_locs, [test] + act_locs, node_index_map) # (n,act,perf,dictionary)
        feas, maxError, numfeas = computeFeas_v1(feeder, [test] + act_locs, A, B, indicMat,substation_name,[test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, load_data, headerpath, modelpath,printCurves) # pass in potential actual loc
        lzn_error_dic[test] = maxError
        markFeas(numfeas, test, graph)
        if feas:
            feas_dic = {}
            feas_dic['act'] = [test]
            feas_dic['perf'] = [test]
            feas_dic['lznErr'] = [lzn_error_dic[test]]
            feas_dic['numfeas'] = [numfeas]
            feas_configs += [feas_dic]        

    heatMapName='heatmap_colocated' + '_' + file_name
    heatMapNames.append(heatMapName)
    nx.nx_pydot.write_dot(graph, heatMapName)
    render('dot', 'png', heatMapName)

    return feas_configs, heatMapNames


def runHeatMapProcess(feeder, set_acts, set_perfs, all_act_locs, perf_nodes, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #compute heatmap (assess feas and lzn error on every node of the feeder, color each red/green on diagram)
    #Then for each act-perf node pair, compute heatmap
    #return list of all "green" configs and the associated lzn errors
    #heatmap shows good places to placed act wrt given perf node, NOT good places to place colocated act-perf node
    #all_act_locs and perf_nodes = lists of node names as strings
    a = 0
    graph = feeder.network
    cur_act_locs = set_acts
    cur_perf_nodes = set_perfs
    heatMapNames = [] # collect heat map names as list of strings
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map, depths)
    lzn_error_run_sum = 0
    feas_configs = []
    
    while a < len(all_act_locs): #outer loop, a = number of actuators
        for act in cur_act_locs: 
            markActLoc(graph, act)
            
        lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
        test_nodes = []
        graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider co-located at substation nodes, node 650 and 651
        
        for node in graphNodes_nosub:
            if node not in cur_act_locs:
                test_nodes.append(node)
        
        for test in test_nodes: #inner loop
            # heatmap color indicates good places to place actuator given chosen loc of perf node (not necessarily colocated)
            indicMat = updateStateSpace(feeder, n, [test] + cur_act_locs, [perf_nodes[a]] + cur_perf_nodes, node_index_map)
            print('evaluating actuator node at ', [test] + cur_act_locs,',\n performance node at ', [perf_nodes[a]] + cur_perf_nodes)
          
            feas, maxError, numfeas = computeFeas_v1(feeder, [test] + cur_act_locs, A, B, indicMat, substation_name,[perf_nodes[a]] + cur_perf_nodes, depths, node_index_map, Vbase_ll, Sbase, load_data, headerpath, modelpath, False)
            lzn_error_dic[test] = maxError
            markFeas(numfeas, test, graph)
            
            if feas:
                feas_dic = {}
                feas_dic['act'] = [test] + cur_act_locs
                feas_dic['perf'] = [test] + cur_perf_nodes
                feas_dic['lznErr'] = [lzn_error_dic[test]]
                feas_dic['numfeas']=[numfeas]
                feas_configs += [feas_dic]        
        
        graph.nodes[perf_nodes[a]]['shape'] = 'square'
        # after generate data for heatmap..
        heatMapName = 'NPP_heatmap_step' + str(a) + '_' + file_name
        heatMapNames.append(heatMapName)
        nx.nx_pydot.write_dot(graph, heatMapName)
        render('dot', 'png', heatMapName)
        a += 1 # place actuator
        
        if a <= len(all_act_locs): # choose actuator and finalize assoc perf node
            cur_act_locs = all_act_locs[0:a] # populate cur_act_locs with subset of all_act_locs
            cur_perf_nodes = perf_nodes[0:a] 
            lzn_error_run_sum += lzn_error_dic[cur_act_locs[-1]][0]
            print('The total max linearization error after '+ str(a) +' actuators have been placed = '+ str(lzn_error_run_sum))
            
        # end of while loop
    return feas_configs, lzn_error_run_sum, heatMapNames


def placeMaxColocActs_stopAtInfeas(feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #place colocated actuators until an infeasible loc is tested, then call find_good_colocated and return 
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map, depths)
    test_nodes = []
    act_locs = []
    printCurves=False # your choice on whether to print PVcurves
    graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
        
    for node in graphNodes_nosub:
        if node not in act_locs:
            test_nodes.append(node)
            
    random.seed(3)  # initialize random num generator so results are reproducable
    while test_nodes:       
        rand_test = random.choice(test_nodes)
        print('evaluating actuator and performance node colocated at ',[rand_test] + act_locs) 
        indicMat = updateStateSpace(feeder, n, [rand_test] + act_locs, [rand_test] + act_locs, node_index_map)
        feas, maxError, numfeas = computeFeas_v1(feeder, [rand_test] + act_locs, A, B, indicMat, substation_name, [rand_test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, load_data, headerpath, modelpath,printCurves)

        if feas:
            act_locs += [rand_test]
            test_nodes.remove(rand_test)
        else:
            feas_configs, heatMapNames = find_good_colocated(feeder, act_locs, node_index_map, substation_name, depths,file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath) # makes a heatmap

            # Believe this img below is redundant as find_good_colocated makes a heatmap already
            #nx.nx_pydot.write_dot(graph, 'CPP_heatmap'+ '_' + file_name)
            #render('dot', 'png', 'CPP_heatmap'+ '_' + file_name)
            break
    
    ff.clear_graph(feeder)
    vis.markActuatorConfig(act_locs, feeder, 'nonauto-CPP')
    return act_locs


def place_max_coloc_acts(seedkey,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #place maximum number of colocated actuators
    #if infeas loc tested, randomly select another test node and continue function run
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(n, feeder, node_index_map, depths)
    test_nodes = []
    act_locs = []
    printCurves = False # your choice on whether to print PVcurves
    graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
        
    for node in graphNodes_nosub:
        if node not in act_locs:
            test_nodes.append(node)
    
    random.seed(seedkey)  # random num generator seed so results are reproducable
    while test_nodes:       
        rand_test = random.choice(test_nodes)
        print('evaluating actuator and performance node colocated at ',[rand_test] + act_locs) 
        indicMat = updateStateSpace(feeder, n, [rand_test] + act_locs, [rand_test] + act_locs, node_index_map)
        feas, maxError, numfeas = computeFeas_v1(feeder, [rand_test] + act_locs, A, B, indicMat, substation_name, [rand_test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, load_data, headerpath, modelpath,printCurves)

        if feas:
            act_locs += [rand_test]
            test_nodes = []
            for node in graphNodes_nosub:
                if node not in act_locs:
                    test_nodes.append(node)
        else:
            test_nodes.remove(rand_test)
    print('Randomly placed ',len(act_locs),' actuators, among n=',len(graphNodes_nosub),' nodes on feeder')
    ff.clear_graph(feeder)
    vis.markActuatorConfig(act_locs, feeder, 'auto-CPP_seed'+str(seedkey))
    return act_locs