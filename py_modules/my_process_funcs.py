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
import my_heatmapSetup_funcs as hm


def eval_config(parmObj,feeder, all_act_locs, perf_nodes, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase):
    #returns whether controllability can be achieved for a given actuator performance node configuration
    #also returns the linearization error associated with the feasibility calculation
    #all_act_locs and perf_nodes = lists of node names as strings
    printCurves = True # your choice on whether to print PVcurves
    graph = feeder.network
    A, B,n = hm.setupStateSpace(parmObj,feeder,depths,file_name)
    assert A.shape==(6*n,6*n), "issue: A (from setupStateSpace) or n have incorrect dims"
    indicMat,indicMat_table,phase_loop_check = hm.updateStateSpace(parmObj,feeder, n, all_act_locs, perf_nodes,file_name)
    if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
        feas, maxError, numfeas,bestF,indicMat = hm.computeFeas_v1(parmObj,feeder, all_act_locs, A, B, indicMat, indicMat_table, substation_name, perf_nodes, depths, node_index_map, Vbase_ll, Sbase,  printCurves,file_name)
    else:
        feas=False
        maxError=1
        numfeas=0
    vis.markActuatorConfig(all_act_locs, feeder, file_name) # create diagram with actuator locs marked
    
    print('Actuator configuration is feasible') if feas else print('Actuator configuration is not feasible')
    return feas, maxError, numfeas,bestF, indicMat

    
def find_good_colocated(parmObj,feeder, set_acts, addon_acts, node_index_map,substation_name, depths, file_name, Vbase_ll, Sbase):
    #almost the same as runheatmap process, but only runs once and shows one heatmap indicating which nodes are good to place a co-located act/perf node
    #return list of all "green" configs and the associated lzn errors
    #act_locs == a list of pre-set colocated act/perf locations --> to evaluate an empty network, pass in act_locs == []
    a = 0 # counter for adding new APNPs
    ff.clear_graph(feeder) # clear any previous modifictions made to graph
    graph = feeder.network
    cur_act_locs = set_acts
    heatMapNames = [] # collect heat map names as list of strings
    A, B,n = hm.setupStateSpace(parmObj,feeder,depths,file_name)
    assert A.shape==(6*n,6*n), "issue: A (from setupStateSpace) or n have incorrect dims"
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
        graphNodes_nosub = hm.remove_subst_nodes(feeder, file_name) # dont consider co-located at substation nodes, node 650 and 651
    
        cand_ctrlType=cand_ctrlTypes[a] # for each step all candidate APNPs will be of this type
        parmObj.set_ctrlTypes([cand_ctrlType]+set_ctrlTypes) # to be used in updateStateSpace
        
        for act in cur_act_locs: # mark placed acts in grey
            vis.markActLoc(graph, act)
            
        for node in graphNodes_nosub: # try placing act/perf at all nodes of the network
            if node not in cur_act_locs:
                test_nodes.append(node)

        for test in test_nodes:
            feas=False # default
            print('evaluating act and perf colocated at ',[test] + cur_act_locs) 
            indicMat,indicMat_table,phase_loop_check = hm.updateStateSpace(parmObj,feeder, n, [test] + cur_act_locs, [test] + cur_act_locs,file_name) # (n,act,perf,dictionary)
            if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
                feas, maxError, numfeas,bestF,indicMat = hm.computeFeas_v1(parmObj,feeder, [test] + cur_act_locs, A, B, indicMat, indicMat_table,substation_name,[test] + cur_act_locs, depths, node_index_map,Vbase_ll, Sbase, printCurves,file_name) # pass in potential actual loc
                lzn_error_dic[test] = maxError
            else:
                maxError=1
                numfeas=0

            vis.markFeas(numfeas, test, graph,phase_loop_check)
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


def runHeatMapProcess(parmObj,feeder, set_acts, set_perfs, addon_acts, addon_perfs, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase):
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
    A, B,n = hm.setupStateSpace(parmObj,feeder, depths,file_name)
    assert A.shape==(6*n,6*n), "issue: A (from setupStateSpace) or n have incorrect dims"
    lzn_error_run_sum = 0
    feas_configs = []
    
    while a < len(addon_acts): #outer loop, a = number of actuators
        for act in cur_act_locs: 
            vis.markActLoc(graph, act)
            
        lzn_error_dic = {} #contains maxLznError for each choice of actuator location with node name as key  
        test_nodes = []
        graphNodes_nosub = hm.remove_subst_nodes(feeder, file_name) # dont consider co-located at substation nodes, node 650 and 651
        
        for node in graphNodes_nosub:
            if node not in cur_act_locs:
                test_nodes.append(node)
        
        for test in test_nodes: #inner loop
            feas=False # default
            # heatmap color indicates good places to place actuator given chosen loc of perf node (not necessarily colocated)          
            print('evaluating actuator node at ', [test] + cur_act_locs,',\n performance node at ', [addon_perfs[a]] + cur_perf_nodes)
            indicMat,indicMat_table,phase_loop_check = hm.updateStateSpace(parmObj,feeder, n, [test] + cur_act_locs, [addon_perfs[a]] + cur_perf_nodes,file_name)
            if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
                feas, maxError, numfeas,bestF,indicMat = hm.computeFeas_v1(parmObj,feeder, [test] + cur_act_locs, A, B, indicMat, indicMat_table, substation_name,[addon_perfs[a]] + cur_perf_nodes, depths, node_index_map, Vbase_ll, Sbase, False,file_name) # false for printing PV curves
                lzn_error_dic[test] = maxError
            else:
                maxError=1
                numfeas=0
            
            vis.markFeas(numfeas, test, graph,phase_loop_check)
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
        vis.write_formatted_dot(graph, heatMapName)

        a += 1 # place actuator
        
        if a <= len(addon_acts): # choose actuator and finalize assoc perf node
            cur_act_locs = addon_acts[0:a]+set_acts # populate cur_act_locs with subset of addon_acts
            cur_perf_nodes = addon_perfs[0:a]+set_acts
            lzn_error_run_sum += lzn_error_dic[cur_act_locs[-1]][0]
            #print('The total max linearization error after '+ str(a) +' actuators have been placed = '+ str(lzn_error_run_sum))
            
        # end of while loop
    return feas_configs, lzn_error_run_sum, heatMapNames

def design_config(parmObj,seedkey,numAct,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase):
    # randomly try groups of numAct actuators until find one that works 
    test_nodes = []
    ctrlTypes=[]

    graphNodes_nosub = hm.remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
    
    for node in graphNodes_nosub:
        test_nodes.append(node)

    random.seed(seedkey)  # initialize random num generator so results are reproducable 
    numtry=0
    while numtry<100: # try a max of 100 configs with numAct actuators
        # choose random set of control types
        if parmObj.get_version()==1:
            ctrlTypes_choosefrom=['PBC','PBC']
        else:
            ctrlTypes_choosefrom=['VWC','VVC']    
            
        rand_ctrlType=random.choices(ctrlTypes_choosefrom,k=numAct)
        ctrlTypes=rand_ctrlType # save into list of control types

        # choose random set of collocated APNP locations
        rand_test = random.sample(test_nodes,numAct) # Return a list of unique elements chosen from the population sequence or set. Used for random sampling without replacement.
        parmObj.set_ctrlTypes(ctrlTypes)
        print('control types=',ctrlTypes)
        print('evaluating actuator and performance node colocated at ',rand_test) 
        feas, maxError,numfeas,bestFparm,indicMat=eval_config(parmObj,feeder, rand_test, rand_test, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase)
        
        if feas:
            break # break out of loop
        else:
            numtry+=1
            
    if not feas:
        print('Algo failed: could not find config with ',numAct,' APNPs')
        bestFparm=-1
    act_locs=rand_test
    return parmObj,act_locs,bestFparm,indicMat

def eval_random_configs(parmObj, seedkey,numAct,numEval,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase):
    # randomly try numEval groups of numAct actuators and returns vector indicating which are feas
    test_nodes = []
    ctrlTypes=[]

    graphNodes_nosub = hm.remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
    
    for node in graphNodes_nosub:
        test_nodes.append(node)

    random.seed(seedkey)  # initialize random num generator so results are reproducable 
    if parmObj.get_version()==1:
        ctrlTypes_choosefrom=['PBC','PBC']
    else:
        ctrlTypes_choosefrom=['VWC','VVC']  
           
    feas_vec=[] # hold 1 if feas, 0 if not
    for i in range(1,numEval): # try a max of 100 configs with numAct actuators
        # choose random set of control types  
        rand_ctrlType=random.choices(ctrlTypes_choosefrom,k=numAct)
        ctrlTypes=rand_ctrlType # save into list of control types

        # choose random set of collocated APNP locations
        rand_test = random.sample(test_nodes,numAct) # Return a list of unique elements chosen from the population sequence or set. Used for random sampling without replacement.
        parmObj.set_ctrlTypes(ctrlTypes)
        print('control types=',ctrlTypes)
        print('evaluating actuator and performance node colocated at ',rand_test) 
        feas, maxError,numfeas,bestFparm,indicMat=eval_config(parmObj,feeder, rand_test, rand_test, node_index_map, substation_name, depths, file_name, Vbase_ll, Sbase)
        
        if feas:
            feas_vec.append(1)
        else: 
            feas_vec.append(0)

    return feas_vec


def placeMaxColocActs_stopAtInfeas(parmObj,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase):
    #place colocated actuators until an infeasible loc is tested, then call find_good_colocated and return 
    graph = feeder.network
    A, B,n = hm.setupStateSpace(parmObj,feeder,depths,file_name)
    assert A.shape==(6*n,6*n), "issue: A (from setupStateSpace) or n have incorrect dims"
    test_nodes = []
    act_locs = []
    ctrlTypes=[]
    printCurves=False # your choice on whether to print PVcurves
    graphNodes_nosub = hm.remove_subst_nodes(feeder, file_name) # dont consider substation nodes, node 650 and 651 for 13NF
    if parmObj.get_version()==1:
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
        indicMat,indicMat_table,phase_loop_check = hm.updateStateSpace(parmObj,feeder, n, [rand_test] + act_locs, [rand_test] + act_locs,file_name)
        if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
            feas, maxError, numfeas,bestF,indicMat = hm.computeFeas_v1(parmObj,feeder, [rand_test] + act_locs, A, B, indicMat, indicMat_table,substation_name, [rand_test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, printCurves, file_name)
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
            feas_configs, heatMapNames = find_good_colocated(parmObj,feeder,[],act_locs, node_index_map, substation_name, depths,file_name, Vbase_ll, Sbase) # makes a heatmap, assume set_acts=[]

            
            (parmObj,feeder, set_acts, addon_acts, node_index_map,substation_name, depths, file_name, Vbase_ll, Sbase)
            
            break
    
    ff.clear_graph(feeder)
    vis.markActuatorConfig(act_locs, feeder, 'nonauto-CPP')
    return act_locs, parmObj


def place_max_coloc_acts(parmObj,seedkey,feeder, file_name, node_index_map, depths, substation_name,Vbase_ll, Sbase, cby_cand=-1):
    #place maximum number of colocated actuators
    #if infeas loc tested, randomly select another test node and continue function run
    # cby_cand is set of all nodes that we want to consider placing an actuator; default set is all nodes without the substation node(s)
    
    graph = feeder.network
    A, B,n = hm.setupStateSpace(parmObj,feeder, depths,file_name)
    assert A.shape==(6*n,6*n), "issue: A (from setupStateSpace) or n have incorrect dims"
    test_nodes = []
    act_locs = []
    ctrlTypes=[]
    printCurves = False # your choice on whether to print PVcurves
    if cby_cand==-1:
        cby_cand=hm.remove_subst_nodes(feeder, file_name) #default set is all nodes without the substation node(s)
        
    print('parmObj.get_version=',parmObj.get_version())
    if parmObj.get_version()==1:
        ctrlTypes_choosefrom=['PBC','PBC']
    else:
        ctrlTypes_choosefrom=['VWC','VVC']
    rand_ctrlType=random.choice(ctrlTypes_choosefrom)
    ctrlTypes.append(rand_ctrlType) # save into list of control types

    for node in cby_cand:
        if node not in act_locs:
            test_nodes.append(node)
    
    random.seed(seedkey)  # random num generator seed so results are reproducable
    while test_nodes:      # while non-empty
        rand_test = random.choice(test_nodes) # choose node randomly among test_nodes
        parmObj.set_ctrlTypes(ctrlTypes)
        print('control types=',ctrlTypes)
        
        print('evaluating actuator and performance node colocated at ',[rand_test] + act_locs) 
        indicMat,indicMat_table,phase_loop_check = hm.updateStateSpace(parmObj,feeder, n, [rand_test] + act_locs, [rand_test] + act_locs,file_name)
        if phase_loop_check:  # disallow configs in which the act and perf node phases are not aligned
            feas, maxError, numfeas,bestF,indicMat = hm.computeFeas_v1(parmObj,feeder, [rand_test] + act_locs, A, B,indicMat, indicMat_table,substation_name, [rand_test] + act_locs, depths, node_index_map,Vbase_ll, Sbase, printCurves, file_name)
        
        else:
            feas=False
            maxError=1
            numfeas=0

        if feas:
            act_locs += [rand_test] # store feas act_loc
            test_nodes = []
            rand_ctrlType=random.choice(ctrlTypes_choosefrom) # choose control for next candidate APNP 
            ctrlTypes.append(rand_ctrlType) # save into list of control types
            for node in cby_cand:
                if node not in act_locs:
                    test_nodes.append(node)
        else:
            test_nodes.remove(rand_test)
    print('Randomly placed ',len(act_locs),' actuators, among n=',len(cby_cand),' possible DER nodes of feeder')
    ff.clear_graph(feeder)
    vis.markActuatorConfig(act_locs, feeder, 'auto-CPP_seed'+str(seedkey))
    return act_locs, parmObj