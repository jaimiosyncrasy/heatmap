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


def setupStateSpace(parmObj,feeder, n,node_index_map, depths):
    #initializes state space matrices A and B
    #n = number of nodes in network
    #feeder = initiaized feeder object
    #node_index_map = dictionary of node indices with node names as keys
    R, X = createRXmatrices_3ph(feeder, node_index_map,depths)
    concat_XR=np.concatenate((X, R), axis = 1)
    
    if parmObj.get_version()==1: # PBC       
        A = np.identity(6*n)
        concat_XR_halfs = np.concatenate(((-1/2) * R, (1/2) * X), axis = 1)
        B = np.concatenate((concat_XR, concat_XR_halfs))
    else: # volt-watt and volt-var
        A = np.zeros((3*n,3*n))
        B = concat_XR # (6n*3n) matrix
    
    return A, B


# correct version
def computeFeas_v1(parmObj,feeder, act_locs, A, B, indicMat, substation_name, perf_nodes, depths, node_index_map, Vbase_ll, Sbase, load_data, headerpath, modelpath, printCurves,file_name):
    node_0 = list(feeder.network.successors(substation_name))
    node_1 = list(feeder.network.successors(node_0[0]))
    z12 = imp.get_total_impedance_from_substation(feeder, node_1[0],depths) # 3 phase, not pu
    B12=np.zeros((3,3)) # TEMPORARY, line susceptance, Yshunt=G+jB

    MYfeas,MYfeasFs,MYnumfeas,MYnumTried,MYnumact,MYbestF,MYindicMat = ctrl.detControlMatExistence(parmObj,feeder, act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map,file_name)
    print('num feas=',MYnumfeas)
    print('num tried=',MYnumTried)

    #lzn_err_max, slopes = lzn.detLznRange(feeder, Vbase_ll, Sbase, z12, B12, act_locs, load_data, headerpath, substation_name, modelpath, depths,printCurves) # usually called by computeFeas
    lzn_err_max=[-1, -1, -1, -1] # workaround, for [PV, QV, Pdel,Qdel] lzn errors

    return MYfeas,lzn_err_max,MYnumfeas,MYbestF,MYindicMat

# workaround version
def computeFeas_v2(feeder, act_locs, perf_nodes, A, B, indicMat):
    random_bit = np.random.choice(2,1)
    if random_bit[0] == 0:
        feas = True
    else:
        feas = False
    maxError = 10
    return feas, maxError

class configParms: # used by updateStateSpace to determine whether each act is PBC, VVC, or VWC
    def __init__(self, default = {}): 
        self.ctrlTypeList = default 
        self.version=-1
    def get_ctrlTypes(self): 
        return self.ctrlTypeList 
    def set_ctrlTypes(self, list): 
        self.ctrlTypeList = list
    def get_version(self): 
        return self.version 
    def set_version(self, ver): 
        self.version = ver
        
def updateStateSpace(parmObj,feeder, n, act_locs, perf_nodes, node_index_map):
    #creates (6n*3n) matrix with 1 at (3i+ph)(3j+ph) for volt-watt control, and (3i+3n+ph)(3j+ph) for volt-var control
    #in the above description, ph is the integer representation (a=0, b=1, c=2) of the phase intersection between the actuator and performance nodes
    #if an actuator and performance node have no phases in common, a warning is printed
    #n = number of nodes in network
    #act_locs = list of actuators locations in network (list of strings)
    #perf_nodes = list of performance nodes 
    #node_index_map = dictionary of node indices for indicMat and F matrix
    if parmObj.get_version()==1: # PBC
        indicMat = np.zeros((6*n,6*n))
    else: # volt-watt and volt-var
        indicMat = np.zeros((6*n,3*n))
    
    ctrlTypeList=parmObj.get_ctrlTypes()
    for i in range(len(act_locs)): 
        act = act_locs[i]
        perf = perf_nodes[i]
        #print('act_locs=',act_locs[i])

        ctrlType=ctrlTypeList[i] # need to have 5 control types, 4 for existing and 1 for the test
        if not(ctrlType=='PBC' or ctrlType=='VVC' or ctrlType=='VWC'):
            raise Exception('Actuator node first 3 chars should be PBC, VVC, or VWC')
        
        act_phases = feeder.busdict[act[4:]].phases # [7:] extracts the YYY bus number from 'XXXbus_YYY'
        perf_phases = feeder.busdict[perf[4:]].phases
        act_index = node_index_map[act] # skip first 3 chars, which is ctrlType
        perf_index = node_index_map[perf]
        
        phase_intrsct = [ph for ph in act_phases if ph in perf_phases]
        if phase_intrsct == []: # disallow configs in which the act and perf node phases are not aligned. Results in kgain=0.0001 and thinks it's feasible
            print('WARNING: act_node ' + act + ' can NOT track perf_node ' + perf + ' --> no common phases')
            phase_loop_check=False
            break # if any actuator of the config has phase mismtach (between act and perf), don't evaluate the config
        else:
            phase_loop_check=True
            for i in range(len(phase_intrsct)):
                if phase_intrsct[i] == 'a':
                    phase_intrsct[i] = 0
                elif phase_intrsct[i] == 'b':
                    phase_intrsct[i] = 1
                elif phase_intrsct[i] == 'c':
                    phase_intrsct[i] = 2

        #print('act=',act,', perf=',perf,', ctrlType=',ctrlType,', phaseItrsct=',phase_intrsct)
        if ctrlType=='PBC':
            for ph in phase_intrsct:
                indicMat[(act_index*3)+ph][(perf_index*3)+ph] = 1
                indicMat[(act_index*3)+(3*n)+ph][(perf_index*3)+(3*n)+ph] = 1  
        elif ctrlType=='VVC':
            for ph in phase_intrsct:
                indicMat[(act_index*3)+ph][(perf_index*3)+ph] = 1
        elif ctrlType=='VWC': #volt-watt control
            for ph in phase_intrsct:
                indicMat[(act_index*3)+(3*n)+ph][(perf_index*3)+ph] = 1   
            
    return indicMat,phase_loop_check

    
def markActLoc(graph, act_loc):
    #changes color of nodes with set actuators to gray
    #graph = networkx graph object (feeder.network)
    #act_loc = node name as string where actuator is placed
    graph.nodes[act_loc]['style'] = 'filled'
    graph.nodes[act_loc]['fillcolor'] = 'indigo'
    graph.nodes[act_loc]['shape'] = 'circle'
    return
        

def markFeas(numfeas, test_act_loc, graph,phase_loop_check):
    #if controllability can be achieved with actuator at test_act_loc then mark green, if only a few controllable configurations (given by thresh_yellowgreen) exist mark yellow, otherwise mark red
    #feas = True or False
    #test_act_loc = node name as string
    #graph = networkx graph object (feeder.network)
    graph.nodes[test_act_loc]['style'] = 'filled'
    graph.nodes[test_act_loc]['shape'] = 'circle'

    thresh_yellowgreen = 15 # you choose
    # for more colors see https://graphviz.org/doc/info/attrs.html
    if not phase_loop_check:
        graph.nodes[test_act_loc]['fillcolor'] = 'firebrick'
    elif numfeas >= thresh_yellowgreen:
        graph.nodes[test_act_loc]['fillcolor'] = 'turquoise'
    elif numfeas >= 1:
        graph.nodes[test_act_loc]['fillcolor'] = 'yellow'
    else:
        graph.nodes[test_act_loc]['fillcolor'] = 'red'
    return


def eval_config(parmObj,feeder, all_act_locs, perf_nodes, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #returns whether controllability can be achieved for a given actuator performance node configuration
    #also returns the linearization error associated with the feasibility calculation
    #all_act_locs and perf_nodes = lists of node names as strings
    printCurves = True # your choice on whether to print PVcurves
    
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(parmObj,feeder,n,node_index_map,depths)
    indicMat,phase_loop_check = updateStateSpace(parmObj,feeder, n, all_act_locs, perf_nodes, node_index_map)
    if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
        feas, maxError, numfeas,bestF,indicMat = computeFeas_v1(parmObj,feeder, all_act_locs, A, B, indicMat, substation_name, perf_nodes, depths, node_index_map, Vbase_ll, Sbase, load_data, headerpath, modelpath, printCurves,file_name)
    else:
        feas=False
        maxError=1
        numfeas=0
    vis.markActuatorConfig(all_act_locs, feeder, file_name) # create diagram with actuator locs marked
    
    print('Actuator configuration is feasible') if feas else print('Actuator configuration is not feasible')
    return feas, maxError, numfeas,bestF, indicMat


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

    
def find_good_colocated(parmObj,feeder, set_acts, addon_acts, node_index_map,substation_name, depths, file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #almost the same as runheatmap process, but only runs once and shows one heatmap indicating which nodes are good to place a co-located act/perf node
    #return list of all "green" configs and the associated lzn errors
    #act_locs == a list of pre-set colocated act/perf locations --> to evaluate an empty network, pass in act_locs == []
    a = 0 # counter for adding new APNPs
    ff.clear_graph(feeder) # clear any previous modifictions made to graph
    graph = feeder.network
    cur_act_locs = set_acts
    heatMapNames = [] # collect heat map names as list of strings
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(parmObj,feeder,n,node_index_map, depths)
    lzn_error_run_sum = 0
    feas_configs = [] 
    printCurves=False # your choice on whether to print PVcurves

    all_ctrlTypes=parmObj.get_ctrlTypes() # format is ctrl types of [set_acts addon_acts]
    set_ctrlTypes=all_ctrlTypes[:len(set_acts)] # get control types of set_acts
    cand_ctrlTypes=all_ctrlTypes[len(set_acts):] # get control types of addon_acts
    #print('set control types=',set_ctrlTypes)
    #print('cand control types=',cand_ctrlTypes)

    while a < len(addon_acts): #outer loop, a = number of actuators to place        
        lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
        test_nodes = []
        graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider co-located at substation nodes, node 650 and 651
    
        cand_ctrlType=cand_ctrlTypes[a] # for each step all candidate APNPs will be of this type
        parmObj.set_ctrlTypes([cand_ctrlType]+set_ctrlTypes) # to be used in updateStateSpace
        
        for act in cur_act_locs: # mark placed acts in grey
            markActLoc(graph, act)
            
        for node in graphNodes_nosub: # try placing act/perf at all nodes of the network
            if node not in cur_act_locs:
                test_nodes.append(node)

        for test in test_nodes:
            feas=False # default
            print('evaluating act and perf colocated at ',[test] + cur_act_locs) 
            indicMat,phase_loop_check = updateStateSpace(parmObj,feeder, n, [test] + cur_act_locs, [test] + cur_act_locs, node_index_map) # (n,act,perf,dictionary)
            if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
                feas, maxError, numfeas,bestF,indicMat = computeFeas_v1(parmObj,feeder, [test] + cur_act_locs, A, B, indicMat,substation_name,[test] + cur_act_locs, depths, node_index_map,Vbase_ll, Sbase, load_data, headerpath, modelpath,printCurves,file_name) # pass in potential actual loc
                lzn_error_dic[test] = maxError
            else:
                maxError=1
                numfeas=0

            markFeas(numfeas, test, graph,phase_loop_check)
            if feas:
                feas_dic = {}
                feas_dic['act'] = [test]+cur_act_locs
                feas_dic['perf'] = [test]+cur_act_locs
                feas_dic['lznErr'] = [lzn_error_dic[test]]
                feas_dic['numfeas'] = [numfeas]
                feas_configs += [feas_dic]

        heatMapName='CPP_heatmap_step' + str(a+1) + '_' + file_name
        heatMapNames.append(heatMapName)
        vis.write_formatted_dot(graph, heatMapName)

        a += 1 # place actuator
        
        if a <= len(addon_acts): # choose actuator and finalize assoc perf node
            cur_act_locs = addon_acts[0:a]+set_acts # populate cur_act_locs with subset of all_act_locs                
            lzn_error_run_sum += lzn_error_dic[cur_act_locs[-1]][0]         
            set_ctrlTypes=[cand_ctrlType]+set_ctrlTypes # update control types for next step
            
    return feas_configs, lzn_error_run_sum, heatMapNames


def runHeatMapProcess(parmObj,feeder, set_acts, set_perfs, addon_acts, addon_perfs, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #compute heatmap (assess feas and lzn error on every node of the feeder, color each red/green on diagram)
    #Then for each act-perf node pair, compute heatmap
    #return list of all "green" configs and the associated lzn errors
    #heatmap shows good places to placed act wrt given perf node, NOT good places to place colocated act-perf node
    #addon_acts and addon_perfs = lists of node names as strings
    a = 0
    graph = feeder.network
    cur_act_locs = set_acts
    cur_perf_nodes = set_perfs
    heatMapNames = [] # collect heat map names as list of strings
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(parmObj,feeder,n, node_index_map, depths)
    lzn_error_run_sum = 0
    feas_configs = []
    
    while a < len(addon_acts): #outer loop, a = number of actuators
        for act in cur_act_locs: 
            markActLoc(graph, act)
            
        lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
        test_nodes = []
        graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider co-located at substation nodes, node 650 and 651
        
        for node in graphNodes_nosub:
            if node not in cur_act_locs:
                test_nodes.append(node)
        
        for test in test_nodes: #inner loop
            feas=False # default
            # heatmap color indicates good places to place actuator given chosen loc of perf node (not necessarily colocated)          
            print('evaluating actuator node at ', [test] + cur_act_locs,',\n performance node at ', [addon_perfs[a]] + cur_perf_nodes)
            indicMat,phase_loop_check = updateStateSpace(parmObj,feeder, n, [test] + cur_act_locs, [addon_perfs[a]] + cur_perf_nodes, node_index_map)
            if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
                feas, maxError, numfeas,bestF,indicMat = computeFeas_v1(parmObj,feeder, [test] + cur_act_locs, A, B, indicMat, substation_name,[addon_perfs[a]] + cur_perf_nodes, depths, node_index_map, Vbase_ll, Sbase, load_data, headerpath, modelpath, False,file_name) # false for printing PV curves
                lzn_error_dic[test] = maxError
            else:
                maxError=1
                numfeas=0
            
            markFeas(numfeas, test, graph,phase_loop_check)
            if feas:
                feas_dic = {}
                feas_dic['act'] = [test] + cur_act_locs
                feas_dic['perf'] = [test] + cur_perf_nodes
                feas_dic['lznErr'] = [lzn_error_dic[test]]
                feas_dic['numfeas']=[numfeas]
                feas_configs += [feas_dic]        
        
        graph.nodes[addon_perfs[a]]['shape'] = 'square'
        # after generate data for heatmap..
        heatMapName = 'NPP_heatmap_step' + str(a+1) + '_' + file_name
        heatMapNames.append(heatMapName)
        write_formatted_dot(graph, heatMapName)

        a += 1 # place actuator
        
        if a <= len(addon_acts): # choose actuator and finalize assoc perf node
            cur_act_locs = addon_acts[0:a]+set_acts # populate cur_act_locs with subset of addon_acts
            cur_perf_nodes = addon_perfs[0:a]+set_acts
            lzn_error_run_sum += lzn_error_dic[cur_act_locs[-1]][0]
            #print('The total max linearization error after '+ str(a) +' actuators have been placed = '+ str(lzn_error_run_sum))
            
        # end of while loop
    return feas_configs, lzn_error_run_sum, heatMapNames


def placeMaxColocActs_stopAtInfeas(parmObj,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #place colocated actuators until an infeasible loc is tested, then call find_good_colocated and return 
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(parmObj,feeder,n, node_index_map, depths)
    test_nodes = []
    act_locs = []
    ctrlTypes=[]
    printCurves=False # your choice on whether to print PVcurves
    graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
    if parmObj.get_version==1:
        ctrlTypes_choosefrom=['PBC','PBC']
    else:
        ctrlTypes_choosefrom=['VWC','VVC']
    rand_ctrlType=random.choice(ctrlTypes_choosefrom)
    ctrlTypes.append(rand_ctrlType) # save into list of control types
    
    for node in graphNodes_nosub:
        if node not in act_locs:
            test_nodes.append(node)
            
    random.seed(3)  # initialize random num generator so results are reproducable
    while test_nodes:       
        rand_test = random.choice(test_nodes)
        parmObj.set_ctrlTypes(ctrlTypes)
        print('control types=',ctrlTypes)
        print('evaluating actuator and performance node colocated at ',[rand_test] + act_locs) 
        indicMat,phase_loop_check = updateStateSpace(parmObj,feeder, n, [rand_test] + act_locs, [rand_test] + act_locs, node_index_map)
        if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
            feas, maxError, numfeas,bestF,indicMat = computeFeas_v1(parmObj,feeder, [rand_test] + act_locs, A, B, indicMat, substation_name, [rand_test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, load_data, headerpath, modelpath,printCurves, file_name)
        else:
            feas=False
            maxError=1
            numfeas=0

        if feas:
            act_locs += [rand_test]
            test_nodes.remove(rand_test)
            rand_ctrlType=random.choice(ctrlTypes_choosefrom)  # choose control for next candidate APNP 
            ctrlTypes.append(rand_ctrlType) # save into list of control types
        else:
            print('Random choice of co-located APNP yields unstable  configuration. Generating heatmap by checking all remaining feeder nodes...')
            feas_configs, heatMapNames = find_good_colocated(parmObj,feeder, act_locs, node_index_map, substation_name, depths,file_name, Vbase_ll, Sbase, load_data, headerpath, modelpath) # makes a heatmap

            break
    
    ff.clear_graph(feeder)
    vis.markActuatorConfig(act_locs, feeder, 'nonauto-CPP')
    return act_locs


def place_max_coloc_acts(parmObj,seedkey,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase, load_data, headerpath, modelpath):
    #place maximum number of colocated actuators
    #if infeas loc tested, randomly select another test node and continue function run
    graph = feeder.network
    n = len(graph.nodes) #number of nodes in network
    A, B = setupStateSpace(parmObj,feeder,n, node_index_map, depths)
    test_nodes = []
    act_locs = []
    ctrlTypes=[]
    printCurves = False # your choice on whether to print PVcurves
    graphNodes_nosub = remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
    if parmObj.get_version==1:
        ctrlTypes_choosefrom=['PBC','PBC']
    else:
        ctrlTypes_choosefrom=['VWC','VVC']
    rand_ctrlType=random.choice(ctrlTypes_choosefrom)
    ctrlTypes.append(rand_ctrlType) # save into list of control types

    for node in graphNodes_nosub:
        if node not in act_locs:
            test_nodes.append(node)
    
    random.seed(seedkey)  # random num generator seed so results are reproducable
    while test_nodes:       
        rand_test = random.choice(test_nodes)
        parmObj.set_ctrlTypes(ctrlTypes)
        print('control types=',ctrlTypes)
        
        print('evaluating actuator and performance node colocated at ',[rand_test] + act_locs) 
        indicMat,phase_loop_check = updateStateSpace(parmObj,feeder, n, [rand_test] + act_locs, [rand_test] + act_locs, node_index_map)
        if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
            feas, maxError, numfeas,bestF,indicMat = computeFeas_v1(parmObj,feeder, [rand_test] + act_locs, A, B, indicMat, substation_name, [rand_test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, load_data, headerpath, modelpath,printCurves, file_name)
        
        else:
            feas=False
            maxError=1
            numfeas=0

        if feas:
            act_locs += [rand_test]
            test_nodes = []
            rand_ctrlType=random.choice(ctrlTypes_choosefrom) # choose control for next candidate APNP 
            ctrlTypes.append(rand_ctrlType) # save into list of control types
            for node in graphNodes_nosub:
                if node not in act_locs:
                    test_nodes.append(node)
        else:
            test_nodes.remove(rand_test)
    print('Randomly placed ',len(act_locs),' actuators, among n=',len(graphNodes_nosub),' nodes on feeder')
    ff.clear_graph(feeder)
    vis.markActuatorConfig(act_locs, feeder, 'auto-CPP_seed'+str(seedkey))
    return act_locs