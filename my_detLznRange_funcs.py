import importlib
import setup_nx # your own module, setup.nx.py
import numpy as np
import math as m
import statistics as st
import cmath
import matplotlib.pyplot as plt 
import itertools
import xlrd
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
import my_heatmapSetup_funcs as hm



def segment_network(feeder, substation_name):
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
            branches += vis.assign_network_branches1(feeder, child)
            
    else:
        branches += [branch_builder]
        return branches
    
    return branches 


def node_index_map_without_substation(graph):
    # make 'special' node index map
    node_index_map = {} 
    t = 0 #initializing first index for node_index_map
    nodes = list(graph.nodes)
    nodes.remove(substation_name)
    
    for node in nodes: #populating node_idex_map
        node_index_map[node] = t
        t += 1
    
    return node_index_map


def trace_error_to_edge(node_name, branch_lst, graph):
    edges_to_update = []
    child_nodes = list(graph.successors(node_name))
    if child_nodes == []:
        edges_to_update += [node_name]
        return edges_to_update
    
    for child in child_nodes:
        cur_branch = find_branch_in_branch_list(child, branch_lst)
        
        if len(cur_branch) >= 2:
            edge = cur_branch[len(cur_branch) - 1: len(cur_branch) - 2: -1][0]
                
        else:
            edge = cur_branch[0]
        
        edges_to_update += trace_error_to_edge(edge, branch_lst, graph)
    
    return edges_to_update

def retrieve_headers(headerpath):
    # retrieve headers from load data file
    load_data_table = pd.read_csv(headerpath, header = None, low_memory = False) # file used to extract table headers (node names) ==> not used to calculate actual loads
    num_cols = len(load_data_table.columns)
    headers = []  
    
    for col_index in range(num_cols)[1:]:
        col = load_data_table[load_data_table.columns[col_index]].tolist()
        my_val = col[0]
        headers += [my_val]
    return headers


def parse_headers(feeder, n, P_vals, Q_vals, headers, node_index_map, modelpath):
    # parse through load data headers to find node names then compute S
    S = np.zeros((3, n))
    S = S.astype('complex128')
    
    workbook_loc = (modelpath) 
    wb = xlrd.open_workbook(workbook_loc) 
    sheet = wb.sheet_by_name('Bus')
    all_phases = []
        
    for r in list(range(sheet.nrows))[1:]:
        all_phases += [sheet.cell_value(r, 0)]
    
    for i in range(len(P_vals)):
        label = headers[i]
        node_name = label[3:]
        node_name = node_name[len(node_name)-4::-1]
        node_name = node_name[len(node_name)::-1]
        node_name = 'bus_' + node_name
        
        find_phase = int(label[-1]) - 1
        phases_per_node = []
        for p in all_phases:
            cur_node = p[len(p)-3::-1]
            cur_node = cur_node[len(cur_node)::-1]
            if node_name[4:] == cur_node:
                phases_per_node += p[len(p)-1:len(p)-2:-1]
        
        phase = phases_per_node[find_phase]
        
        if phase == 'a':
            phase_index = 0
        elif phase =='b':
            phase_index = 1
        elif phase == 'c':
            phase_index = 2
                
        node_index = node_index_map[node_name]
        S[phase_index][node_index] = complex(P_vals[i], Q_vals[i])/(5000/3)
    return S


def find_edge_nodes(feeder):
    # Find edge nodes
    graph = feeder.network
    edge_nodes = []
    
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            edge_nodes += [node]
    return edge_nodes

#def calc_losses_for_transformer():
    
    
    
def calc_losses_for_line(node_index_map, feeder, edge_node, Iforks, S, Sloss, Sact, V):
    cur_index = node_index_map[cur_node]
    pred_index = node_index_map[pred_node[0]]
            
    V_drop = np.dot(z_3by3, np.array([ia_prev, ib_prev, ic_prev]))
    V[0][cur_index] = V[0][child_index] + V_drop[0]
    V[1][cur_index] = V[1][child_index] + V_drop[1]
    V[2][cur_index] = V[2][child_index] + V_drop[2]

    ia_conj = (S[0][cur_index] + Sact[cur_index]) / (V[0][cur_index])
    ib_conj = (S[1][cur_index] + Sact[cur_index]) / (V[1][cur_index])
    ic_conj = (S[2][cur_index] + Sact[cur_index]) / (V[2][cur_index])
    # take complex conjugate:
    ia_load = complex(ia_conj.real, (-1) * ia_conj.imag)
    ib_load = complex(ib_conj.real, (-1) * ib_conj.imag)
    ic_load = complex(ic_conj.real, (-1) * ic_conj.imag)

    ia, ib, ic = ia_load + ia_prev, ib_load + ib_prev, ic_load + ic_prev 
    impedance = feeder.network.get_edge_data(pred_node[0], cur_node, default=None)['connector']
    z_3by3 = impedance.Z if isinstance(impedance, line) else np.zeros((3,3))
    V_drop = np.dot(z_3by3, np.array([ia, ib, ic]))
    conj_ia = complex(ia.real, (-1) * ia.imag)
    conj_ib = complex(ib.real, (-1) * ib.imag)
    conj_ic = complex(ic.real, (-1) * ic.imag)
    I_3phase_conj = np.array([conj_ia, conj_ib, conj_ic])
    Sloss_3phase = np.dot(I_3phase_conj, V_drop)
    Sloss[pred_index][cur_index] = Sloss_3phase
                
    cur_node = pred_node[0]
    cur_node_children = list(graph.successors(cur_node))
    pred_node = list(graph.predecessors(cur_node))
    child_index = cur_index
    ia_prev, ib_prev, ic_prev = ia, ib, ic


def calc_line_losses_to_fork(node_index_map, feeder, edge_node, Iforks, S, Sloss, Sact, V, Zbase):
    graph = feeder.network
    pred_node = list(graph.predecessors(edge_node))
    cur_index = node_index_map[edge_node]
    pred_index = node_index_map[pred_node[0]]
    
    ia_conj = (S[0][cur_index] + Sact[cur_index]) / V[0][cur_index] # all in pu
    ib_conj = (S[1][cur_index] + Sact[cur_index]) / V[1][cur_index]
    ic_conj = (S[2][cur_index] + Sact[cur_index]) / V[2][cur_index]
    # take complex conjugate:
    ia_load = complex(ia_conj.real, (-1) * ia_conj.imag)
    ib_load = complex(ib_conj.real, (-1) * ib_conj.imag)
    ic_load = complex(ic_conj.real, (-1) * ic_conj.imag)
            
    ia = ia_load + sum([i[0] for i in Iforks[edge_node]])
    ib = ib_load + sum([i[1] for i in Iforks[edge_node]])
    ic = ic_load + sum([i[2] for i in Iforks[edge_node]])
    
    impedance = feeder.network.get_edge_data(pred_node[0], edge_node, default=None)['connector']
    z_3by3 = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3)) # not pu
    z_3by3 = z_3by3/Zbase # convert to pu 
    
    V_drop = np.dot(z_3by3, np.array([ia, ib, ic]))
    conj_ia = complex(ia.real, (-1) * ia.imag)
    conj_ib = complex(ib.real, (-1) * ib.imag)
    conj_ic = complex(ic.real, (-1) * ic.imag)
    I_3phase_conj = np.array([conj_ia, conj_ib, conj_ic])
    Sloss_3phase = np.dot(I_3phase_conj, V_drop)
    Sloss[pred_index][cur_index] = Sloss_3phase
        
    cur_node = pred_node[0]
    cur_node_children = list(graph.successors(cur_node))
    pred_node = list(graph.predecessors(cur_node))
    child_index = node_index_map[cur_node_children[0]]
    ia_prev, ib_prev, ic_prev = ia, ib, ic

    while len(cur_node_children) == 1 and pred_node != []: 
        cur_index = node_index_map[cur_node]
        pred_index = node_index_map[pred_node[0]]
            
        V_drop = np.dot(z_3by3, np.array([ia_prev, ib_prev, ic_prev]))   
        V[0][cur_index] = V[0][child_index] + V_drop[0]
        V[1][cur_index] = V[1][child_index] + V_drop[1]
        V[2][cur_index] = V[2][child_index] + V_drop[2]

        ia_conj = (S[0][cur_index] + Sact[cur_index]) / (V[0][cur_index])
        ib_conj = (S[1][cur_index] + Sact[cur_index]) / (V[1][cur_index])
        ic_conj = (S[2][cur_index] + Sact[cur_index]) / (V[2][cur_index])
        # take complex conjugate:
        ia_load = complex(ia_conj.real, (-1) * ia_conj.imag)
        ib_load = complex(ib_conj.real, (-1) * ib_conj.imag)
        ic_load = complex(ic_conj.real, (-1) * ic_conj.imag)

        ia, ib, ic = ia_load + ia_prev, ib_load + ib_prev, ic_load + ic_prev 
        impedance = feeder.network.get_edge_data(pred_node[0], cur_node, default=None)['connector']
        z_3by3 = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3)) # not pu
        z_3by3 = z_3by3/Zbase # convert to pu 
        V_drop = np.dot(z_3by3, np.array([ia, ib, ic]))
        conj_ia = complex(ia.real, (-1) * ia.imag)
        conj_ib = complex(ib.real, (-1) * ib.imag)
        conj_ic = complex(ic.real, (-1) * ic.imag)
        I_3phase_conj = np.array([conj_ia, conj_ib, conj_ic])
        Sloss_3phase = np.dot(I_3phase_conj, V_drop)
        Sloss[pred_index][cur_index] = Sloss_3phase
                
        cur_node = pred_node[0]
        cur_node_children = list(graph.successors(cur_node))
        pred_node = list(graph.predecessors(cur_node))
        child_index = cur_index
        ia_prev, ib_prev, ic_prev = ia, ib, ic
    
    return S, V, ia_prev, ib_prev, ic_prev, cur_node, z_3by3, child_index


def compute_line_losses_multiphase(feeder, P_vals, Q_vals, act_locs, Sbase, Zbase, headerpath, substation_name, modelpath, depths, lb_mode = True):
    # P_vals and Q_vals not pu
    graph = feeder.network
    n = len(graph.nodes) 
    percent_V_drop = .05
    P = [p/Sbase for p in P_vals] # converting to pu
    Q = [q/Sbase for q in Q_vals] # converting to pu
    branch_list = segment_network(feeder, substation_name)       
    status = 'unsolved'
    run_counter = 0
    node_index_map = hm.createNodeIndexMap(feeder)
    
    # retrieve headers from load data file
    headers = retrieve_headers(headerpath)
    
    # parse through load data headers to find node names
    S = parse_headers(feeder, n, P, Q, headers, node_index_map, modelpath)
    Sact = np.zeros(n)
    
    # set actuators to extreme values for lower or upper bound
    if lb_mode:
        act_cap = -(50/Sbase) #per phase
    else:
        act_cap = (50/Sbase) #per phase
    for act in act_locs:
        index = node_index_map[act]
        Sact[index] = act_cap
        
    # Find edge nodes:
    edge_nodes = find_edge_nodes(feeder)
    
    edge_nodes.remove('bus_634')
    edge_nodes += ['bus_633']
    
    # Find impedances from edge nodes to substation
    z_edges_to_sub = np.zeros((3, len(edge_nodes)), dtype = complex)
    i = 0
    
    for edge in edge_nodes:
        cur_z = imp.get_total_impedance_from_substation(feeder, edge, depths)
        z_edges_to_sub[0][i] = cur_z[0][0] 
        z_edges_to_sub[1][i] = cur_z[1][1] 
        z_edges_to_sub[2][i] = cur_z[2][2]
        i += 1
    
    #zmax = max([max(z_edges_to_sub[0]), max(z_edges_to_sub[1]), max(z_edges_to_sub[2])])
    #Vest = np.empty((3, n), dtype = complex)
    Vest = np.zeros((3, n), dtype = complex) # Vest is proportional to distance from substaion.
    percent_V_drops = percent_V_drop * np.ones((3, n))
    
    while status == 'unsolved' and run_counter < 1:
        loop_status = 'unbroken'
        run_counter += 1
    
        # estimate voltages at edge nodes
        #for i in range(len(edge_nodes)):
            #cur_edge = edge_nodes[i]
            #cur_index = node_index_map[cur_edge]
            #cur_za, cur_zb , cur_zc = z_edges_to_sub[0][i], z_edges_to_sub[1][i], z_edges_to_sub[2][i]
            #va = 1 - (percent_V_drops[0][cur_index] * (cur_za / zmax)) # per unit
            #vb = 1 - (percent_V_drops[1][cur_index] * (cur_zb / zmax))
            #vc = 1 - (percent_V_drops[2][cur_index] * (cur_zc / zmax))
            #da = -(percent_V_drops[0][cur_index] * (cur_za / zmax)) # per unit
            #db = -(percent_V_drops[1][cur_index] * (cur_zb / zmax))
            #dc = -(percent_V_drops[2][cur_index] * (cur_zc / zmax))
            #Vest[0][cur_index] = complex(va * m.cos(da), va * m.sin(da))
            #Vest[1][cur_index] = complex(vb * m.cos(db), vb * m.sin(db))
            #Vest[2][cur_index] = complex(vc * m.cos(dc), vc * m.sin(dc))
            
        Vest[0][3] = complex(.9859 * m.cos((m.pi/180)*(-121.2784)), .9859 * m.sin(-121.278*(m.pi/180)))
        Vest[1][3] = complex( .9682* m.cos(118.0922*(m.pi/180)), .9682 * m.sin(118.0922*(m.pi/180)))
        Vest[2][3] = complex( .1113* m.cos(-2.004*(m.pi/180)), .1113 * m.sin(-2.004*(m.pi/180)))
        Vest[0][5] = complex( .9871* m.cos(-121.247*(m.pi/180)), .9871* m.sin(-121.247*(m.pi/180)))
        Vest[1][5] = complex(.973 * m.cos(118.008*(m.pi/180)), .973 * m.sin(118.008*(m.pi/180)))
        Vest[2][5] = complex( 1* m.cos(-.0006*(m.pi/180)), 1 * m.sin(-.0006*(m.pi/180)))
        Vest[0][11] = complex(.9964 * m.cos(-121.4939*(m.pi/180)),  .9964* m.sin(-121.4939*(m.pi/180)))
        Vest[1][11] = complex( .9401* m.cos(116.636*(m.pi/180)), .9401 * m.sin(116.636*(m.pi/180)))
        Vest[2][11] = complex( .9689* m.cos(-3.69*(m.pi/180)), .9689 * m.sin(-3.69*(m.pi/180)))
        Vest[0][0] = complex( .9393* m.cos(116.54),.9393  * m.sin(116.54))
        Vest[1][0] = complex( .9401* m.cos(116.6367*(m.pi/180)),  .9401* m.sin(116.6367*(m.pi/180)))
        Vest[2][0] = complex(.9819 * m.cos(-1.59*(m.pi/180)), .9819 * m.sin(-1.59*(m.pi/180)))
        Vest[0][8] = complex( .9701* m.cos(-3.69*(m.pi/180)), .9701 * m.sin(-3.69*(m.pi/180)))
        Vest[1][8] = complex( .9401* m.cos(116.6367*(m.pi/180)),  .9401* m.sin(116.6367*(m.pi/180)))
        Vest[2][8] = complex(.9701 * m.cos(-3.69*(m.pi/180)),  .9701* m.sin(-3.69*(m.pi/180)))
        Vest[0][10] = complex(.9984* m.cos(-121.63*(m.pi/180)),.9984  * m.sin(-121.63*(m.pi/180)))
        Vest[1][10] = complex( .9392* m.cos(116.6192*(m.pi/180)), .9392 * m.sin(116.6192*(m.pi/180)))
        Vest[2][10] = complex( .9701* m.cos(-3.69*(m.pi/180)), .9701 * m.sin(-3.69*(m.pi/180)))
        

        V = Vest
        Sloss = np.zeros((n,n), dtype = complex) # will populate 
        Vforks = {}
        Iforks = {}
        active_branches = []
        
        for branch in branch_list:
            for edge in edge_nodes:
                if edge in branch:
                    active_branches += [branch]
                    break
                         
            edge = branch[-1]
            Iforks[edge] = []
            Vforks[edge] = []
        
        for _ in range(len(branch_list)):
            if loop_status == 'broken':
                break
                
            cur_branch = active_branches[0]
            edge = cur_branch[-1]
            # all outputs in pu
            S, V, ia_prev, ib_prev, ic_prev, cur_node, z_3by3, child_index = calc_line_losses_to_fork(node_index_map, feeder, edge, Iforks, S, Sloss, Sact, V, Zbase)
            
            if cur_node == substation_name:
                #end of function run
                Seq = sum(sum(Sloss)) + sum(sum(S)) + sum(Sact)
                #Seq = sum(sum(S)) + sum(Sact)
                Peq = Seq.real
                Qeq = Seq.imag
                print('Number of iterations performed: ' + str(run_counter))
                return Peq, Qeq
            
            V_drop = np.dot(z_3by3, np.array([ia_prev, ib_prev, ic_prev]))
            Vfork = [edge, V[0][child_index] + V_drop[0]]
            Vfork += [V[1][child_index] + V_drop[1]]
            Vfork += [V[2][child_index] + V_drop[2]]
            Vforks[cur_node] += [Vfork]
            Iforks[cur_node] += [[ia_prev, ib_prev, ic_prev]]
            active_branches.remove(cur_branch)
            
            if len(active_branches) == 0:
                for key, val in Vforks.items():
                    if loop_status == 'broken':
                        break
                    
                    # check to see if all child branches of fork have been iterated through
                    if len(val) > 0 and len(list(graph.successors(key))) == len(val):
                        for v_outer in val:
                            if key == 'bus_632':
                                break                                                                                  
                            if loop_status == 'broken':
                                break
                                
                            for v_inner in val:
                                check_a_real = (v_outer[1].real - v_inner[1].real) / v_outer[1].real
                                check_a_imag = (v_outer[1].imag - v_inner[1].imag) / v_outer[1].imag
                                check_b_real = (v_outer[2].real - v_inner[2].real) / v_outer[2].real
                                check_b_imag = (v_outer[2].imag - v_inner[2].imag) / v_outer[2].imag
                                check_c_real = (v_outer[3].real - v_inner[3].real) / v_outer[3].real
                                check_c_imag = (v_outer[3].imag - v_inner[3].imag) / v_outer[3].imag
                                all_checks = [check_a_real, check_a_imag, check_b_real, check_b_imag, check_c_real, check_c_imag]
                                
                                #checkVclose = [abs(check_a.real) < .1 and abs(check_a.imag) < .1] # check < than 10% different
                                #checkVclose += [abs(check_b.real) < .1 and abs(check_b.imag) < .1]
                                #checkVclose += [abs(check_c.real) < .1 and abs(check_c.imag) < .1]
                                
                                # check not more than 10% different
                                large_error_indxs = [i for i in [0, 0, 1, 1, 2, 2] if abs(all_checks[i]) > .1]   
                                for i in large_error_indxs:
                                    if all_checks[i] < 0:
                                        high_v = v_inner
                                    else:
                                        high_v = v_outer                                    
                                    lower_edge_v = trace_error_to_edge(high_v[0], branch_list, graph)
                                    for edge in lower_edge_v:
                                        cur_index = node_index_map[edge]
                                        #print('w')
                                        percent_V_drops[i][cur_index] -= .5
                                    
                                if len(large_error_indxs) > 0:
                                    #print(Vforks['bus_632'])
                                    
                                    #for n in graph.nodes:
                                        #i = node_index_map[n]
                                        #v = [V[0][i]]
                                        #v += [V[1][i]]
                                        #v += [V[2][i]]
                                    print(key)
                                        #print(v)
                        
                                    loop_status = 'broken'
                                    break                            
                            
                            #for v_inner in val:
                                #checkVclose = [(v_outer[0] - v_inner[0]) / v_outer[0]] # check not more than 10% different
                                #checkVclose += [(v_outer[1] - v_inner[1]) / v_outer[1]]
                                #checkVclose += [(v_outer[2] - v_inner[2]) / v_outer[2]]
                                #print(key + ": " + str(checkVclose))
                        if key == 'bus_632':
                            key_index = node_index_map[key]
                            v = [v for v in val if v[0] == 'bus_671'][0]
                            V[0][key_index] = v[1]
                            V[1][key_index] = v[2]
                            V[2][key_index] = v[3]                            
                            for branch in branch_list:
                                if key in branch:
                                    active_branches += [branch]
                                    Vforks[key] = []
                                    break

                        if loop_status == 'unbroken' and key != 'bus_632':
                            key_index = node_index_map[key]
                            V[0][key_index] = np.mean([v[1] for v in val])
                            V[1][key_index] = np.mean([v[2] for v in val])
                            V[2][key_index] = np.mean([v[3] for v in val])                           
                            for branch in branch_list:
                                if key in branch:
                                    active_branches += [branch]
                                    Vforks[key] = []
                                    break


#--------------------------------------------------------------------------------------
        
def find_timestep(PQ_extreme, run_sum, Qstart_index, data_table):
    time_col = data_table[data_table.columns[0]].tolist()
    num_cols = len(data_table.columns)
    
    for i in range(len(run_sum)):
        
        if PQ_extreme == run_sum[i]:
            time = time_col[i]
            extreme_row = [] 
            
            for col_index in range(num_cols)[1:]:
                col = data_table[data_table.columns[col_index]].tolist()
                my_val = col[i]
                extreme_row += [my_val]
            break
    
    extreme_row_P = extreme_row[0:Qstart_index - 1]
    extreme_row_Q = extreme_row[Qstart_index - 1:]
    return time, extreme_row_P, extreme_row_Q


def computePQsweep_timesteps(feeder, load_data):
    # load_data = name as string of csv file with time varying load data for a network 
    data_table = pd.read_csv(load_data)
    num_cols = len(data_table.columns)
    time_col = data_table[data_table.columns[0]].tolist()
    Qstart_index = round(((num_cols - 1) / 2) + 1)
    run_sum_P = data_table[data_table.columns[1]].tolist()
    run_sum_Q = data_table[data_table.columns[Qstart_index]].tolist()
    
    for i in range(num_cols)[2:Qstart_index]:
        col = data_table[data_table.columns[i]].tolist()
        rum_sum_P = list(map(add, col, run_sum_P))
        
    for i in range(num_cols)[(Qstart_index + 1):]:
        col = data_table[data_table.columns[i]].tolist()
        rum_sum_Q = list(map(add, col, run_sum_Q))
    
    #max/min over timesteps
    Psweep_lb = min(run_sum_P)
    Psweep_ub = max(run_sum_P)
    Qsweep_lb = min(run_sum_Q)
    Qsweep_ub = max(run_sum_Q)
    
    print('P_lb(kW) = '+str(Psweep_lb))
    print('P_ub(kW) = '+str(Psweep_ub))
    print('Q_lb(kVAR) = '+str(Qsweep_lb))
    print('Q_ub(kVAR) = '+str(Qsweep_ub))
    
    time_P_lb, P_lb_rowP, P_lb_rowQ = find_timestep(Psweep_lb, run_sum_P, Qstart_index, data_table) 
    time_P_ub, P_ub_rowP, P_ub_rowQ = find_timestep(Psweep_ub, run_sum_P, Qstart_index, data_table)
    time_Q_lb, Q_lb_rowP, Q_lb_rowQ = find_timestep(Qsweep_lb, run_sum_Q, Qstart_index, data_table)
    time_Q_ub, Q_ub_rowP, Q_ub_rowQ = find_timestep(Qsweep_ub, run_sum_Q, Qstart_index, data_table)
    
    P_lb_results = [time_P_lb, P_lb_rowP, P_lb_rowQ]
    P_ub_results = [time_P_ub, P_ub_rowP, P_ub_rowQ]
    Q_lb_results = [time_Q_lb, Q_lb_rowP, Q_lb_rowQ]
    Q_ub_results = [time_Q_ub, Q_ub_rowP, Q_ub_rowQ]
    
    print('Timesteps for Extreme Load Values:')
    print('P_lb = ' + str(time_P_lb))
    print('P_ub = ' + str(time_P_ub))
    print('Q_lb = ' + str(time_Q_lb))
    print('Q_ub = ' + str(time_Q_ub))
    
    return P_lb_results, P_ub_results, Q_lb_results, Q_ub_results


def computePQsweep_losses(feeder, act_locs, Sbase, Zbase, P_lb_results, P_ub_results, Q_lb_results, Q_ub_results, headerpath, substation_name, modelpath, depths):
    time_P_lb, P_lb_rowP, P_lb_rowQ = P_lb_results[0], P_lb_results[1], P_lb_results[2]
    time_P_ub, P_ub_rowP, P_ub_rowQ = P_ub_results[0], P_ub_results[1], P_ub_results[2]
    time_Q_lb, Q_lb_rowP, Q_lb_rowQ = Q_lb_results[0], Q_lb_results[1], Q_lb_results[2]
    time_Q_ub, Q_ub_rowP, Q_ub_rowQ = Q_ub_results[0], Q_ub_results[1], Q_ub_results[2]
    
    # accounting for line losses:
    Psweep_lb, _ = compute_line_losses_multiphase(feeder, P_lb_rowP, P_lb_rowQ, act_locs, Sbase, Zbase, headerpath, substation_name, modelpath, depths, lb_mode = True)
    Psweep_ub, _ = compute_line_losses_multiphase(feeder, P_ub_rowP, P_ub_rowQ, act_locs, Sbase, Zbase, headerpath, substation_name, modelpath, depths, lb_mode = False)
    _, Qsweep_lb = compute_line_losses_multiphase(feeder, Q_lb_rowP, Q_lb_rowQ, act_locs, Sbase, Zbase, headerpath, substation_name, modelpath, depths, lb_mode = True)
    _, Qsweep_ub = compute_line_losses_multiphase(feeder, Q_ub_rowP, Q_ub_rowQ, act_locs, Sbase, Zbase, headerpath, substation_name, modelpath, depths, lb_mode = False)
    
    PQ_bounds = [Psweep_lb, Psweep_ub, Qsweep_lb, Qsweep_ub]
    return PQ_bounds


# Solve fwd-bwd sweep single phase
def solveFwdBwdSweep_2bus(R12, X12, V1, P2, Q2):
    # Initialization
    #print('~~~~~~~ Starting FBS Method for Solving PF')

    # Givens: z12, V1, P2, Q2
    S2 = complex(P2, Q2)  # per unit
    z12 = complex(R12, X12)  # per unit
    Vs = V1  # per unit

    # Init Cond
    V1 = []
    V2 = []
    Vconv = []
    Vnom = Vs  # to check convergence

    tol = 0.0001
    k = 0
    V1.append(0)
    V2.append(0)
    Vconv.append([0, 0])

    '''First Iteration'''
    k += 1

    # Fwd Sweep
    V1.append(Vs)
    V2.append(Vs)

    # Check convergence:
    Vconv.append([abs((abs(V1[k]) - abs(V1[k - 1]))) / Vnom, \
                  abs((abs(V2[k]) - abs(V2[k - 1]))) / Vnom])
    # Backward sweep
    I12 = np.conj(S2 / V2[k])

    '''Iterative Part'''
    while any(node >= tol for node in Vconv[k]):  # break when all nodes less than tol
        k += 1  # new iteration
        # Fwd sweep
        V1.append(V1[k - 1])  # same as prev iter ZERO?
        V2.append(Vs - (z12 * I12))

        # Check convergence:
        Vconv.append([abs((abs(V1[k]) - abs(V1[k - 1]))) / Vnom, \
                      abs((abs(V2[k]) - abs(V2[k - 1]))) / Vnom])

        # print(Vconv) uncomment when debugging

        # Backward sweep
        I12 = np.conj(S2 / V2[k])

        if len(Vconv) > 30:
            print('Didnt converge')
            break  # break out of loop of iter too many times

    '''Output Results'''
    #print('~~~~~~~ PF Results: ')
    Vsoln = [V1[-1], V2[-1]]  # /Vs to put into pu
    #print(Vsoln)
    convergedIfZero = Vconv[-1]
    #print(convergedIfZero)
    numIter = len(Vconv) - 1  # -1 because Vconv initialized at zero
    #print(numIter)
    #print('~~~~~~~ Finished FBS Method for SOlving PF')

    '''Polar to rect conversion for testing/probing'''
    mag = [abs(ele) for ele in Vsoln]
    ang = [np.degrees(cmath.phase(ele)) for ele in Vsoln]
    Vsoln_polarCoor = [mag, ang]  # Vpu, deg
    #print('mag:',Vsoln_polarCoor[0])
    #print('ang:', Vsoln_polarCoor[1])
    V2 = abs(Vsoln[1])
    del2 = np.degrees(cmath.phase(Vsoln[1]))

    return V2, del2 # end of solve PF on 2-bus

# Solve fwd-bwd sweep 3-phase
def solveFwdBwdSweep_2bus_3ph(R12, X12, B12,Vs, P2, Q2):
    # Initialization
    #print('~~~~~~~ Starting FBS Method for Solving PF')

    # Givens: z12, V1, P2, Q2
    S2=P2+Q2*1j # per unit
    z12 = R12+X12*1j  # per unit
    Vs=np.transpose(Vs) # make 1x3 vector, Vs is per unit

    # Init Cond
    V1,V2,Vconv=(np.empty((0,3)) for i in range(3))
    Vnom = Vs  # to check convergence
    
    tol = 0.0001
    k = 0
    V1=np.append(V1,np.zeros((1,3)),axis=0)
    V2=np.append(V2,np.zeros((1,3)),axis=0)
    Vconv=np.append(Vconv,np.zeros((2,3)),axis=0)

   #  '''First Iteration'''
    k += 1
    
    # Fwd Sweep
    V1=np.append(V1,Vs,axis=0) # concatenate along rows
    V2=np.append(V2,Vs,axis=0)

     # Check convergence:
    a=np.array(np.divide(abs((abs(V1[k]) - abs(V1[k - 1]))),Vnom))
    b=np.array(np.divide(abs((abs(V2[k]) - abs(V2[k - 1]))),Vnom))
    c=np.concatenate((a, b), axis=0)
    Vconv=np.append(Vconv,c,axis=0)
                   
     # Backward sweep
    I12 = np.conj(np.divide(np.transpose(S2),V2[k]))-np.dot(np.conj(B12*1j),V2[k])
    #print('I12=',I12)
    #print(np.concatenate((Vconv[2*k], Vconv[2*k+1])))
    
    '''Iterative Part'''
    
    while any(node >= tol for node in np.concatenate((Vconv[2*k], Vconv[2*k+1]))):  # break when all nodes less than tol
        k += 1  # new iteration
         # Fwd sweep
        V1=np.append(V1,[V1[k-1]],axis=0)# same as prev iter ZERO?
        V2=np.append(V2,Vs - np.dot(I12,z12),axis=0) # np.dot is matrix mult
        #print('V1=',V1)
        #print('V2=',V2)

        # LEFTOFF here
#         # Check convergence:
        a=np.array(np.divide(abs((abs(V1[k]) - abs(V1[k - 1]))),Vnom))
        b=np.array(np.divide(abs((abs(V2[k]) - abs(V2[k - 1]))),Vnom))
        c=np.concatenate((a, b), axis=0)
        Vconv=np.append(Vconv,c,axis=0)

        #print('Vconv=',Vconv) # uncomment when debugging
    
         # Backward sweep
        I12 = np.conj(np.divide(np.transpose(S2),V2[k]))-np.dot(np.conj(B12*1j),V2[k])

        if len(Vconv) > 30:
            print('Didnt converge')
        break  # break out of loop of iter too many times

    '''Output Results'''
    #print('~~~~~~~ PF Results: ')
    V1soln=V1[-1] #end
    V2soln=V2[-1]
    convergedIfZero = Vconv[-1]
   # print('convergedIfZero=',convergedIfZero)
    numIter = len(Vconv) - 1  # -1 because Vconv initialized at zero
    #print('num iterations=',numIter)
    #print('~~~~~~~ Finished FBS Method for Solving PF')

    '''Polar to rect conversion for testing/probing'''
    mag = [[abs(ele) for ele in V2soln]]
    ang = [[np.degrees(cmath.phase(ele)) for ele in V2soln]]
    Vsoln_polarCoor = [mag, ang]  # Vpu, deg
    #print('V2mag mag (Vpu):',Vsoln_polarCoor[0])
    #print('V2ang ang (degrees):', Vsoln_polarCoor[1])

    V2=np.transpose(mag)
    del2=np.transpose(ang)
    return V2, del2 # Vpu, degrees


def makePVcurve(sweep_lb, sweep_ub, Sbase, Vbase, R12, X12, V1):
    numPts = 20
    P12 = np.linspace(sweep_lb, sweep_ub, numPts)
    Q12pu = m.tan(m.acos(.9))
    Q12 = Q12pu * Sbase
    trueV2 = np.zeros(numPts)
    trueDel2 = np.zeros(numPts)
    lznV2 = np.zeros(numPts)
    lznDel2 = np.zeros(numPts)
    solns = {}
    
    for i in range(len(P12)):
        # create true curve
        a, b = solveFwdBwdSweep_2bus(R12, X12, V1, P12[i], Q12) # inputs are 3x3 matrices
        trueV2[i] = a
        trueDel2[i] = b
        
        # create lzn curve
        V2sq = (V1**2) - (2*R12*P12[i]) - (2*X12*Q12)
        V2 = V2sq**(1/2) # take sqrt
        delta2 = 0 - (((X12*P12[i])-(R12*Q12))/(V1*V2))
        lznV2[i] = V2
        lznDel2[i] = (180/m.pi)*delta2
    
    plt.figure(1)
    plt.plot(P12/Sbase, lznV2/Vbase,'r', label = 'linearization')
    plt.plot(P12/Sbase, trueV2/Vbase,'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('V2, pu')
    plt.title('True P-V Curve and Linearization Curve')
    plt.legend()
    plt.savefig('True_PV_Curve_and_Linearization_Curve.png')
    
    plt.figure(2)
    plt.plot(P12/Sbase, lznDel2,'r', label = 'linearization')
    plt.plot(P12/Sbase, trueDel2,'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('Delta2, degrees')
    plt.title('True P-Del Curve and Linearization Curve')
    plt.legend()
    plt.savefig('True_P-Del_Curve_and_Linearization_Curve.png')
    
    solns['trueV2'] = trueV2
    solns['trueDel2'] = trueDel2 
    solns['lznV2'] = lznV2 
    solns['lznDel2'] = lznDel2
    
    return P12, solns # end of makePVcurve
    
    
def makeQVcurve(Sweep_lb, Sweep_ub, Sbase, Vbase, R12, X12, V1):
    numPts = 20
    Q12 = np.linspace(Sweep_lb, Sweep_ub, numPts)
    P12pu = m.tan(m.acos(0.9))  
    P12 = P12pu * Sbase
    trueV2 = np.zeros(numPts)
    trueDel2 = np.zeros(numPts)
    lznV2 = np.zeros(numPts)
    lznDel2 = np.zeros(numPts)
    solns = {}
    
    for i in range(len(Q12)):
        a, b = solveFwdBwdSweep_2bus(R12, X12, V1, P12, Q12[i])
        trueV2[i] = a
        trueDel2[i] = b
        V2sq = (V1**2) - (2*R12*P12) - (2*X12*Q12[i])
        V2 = V2sq**(1/2)
        delta2 = 0 - (((X12*P12)-(R12*Q12[i]))/(V1*V2))
        lznV2[i] = V2
        lznDel2[i] = (180/m.pi)*delta2
        
    plt.figure(3)
    plt.plot(Q12/Sbase, lznV2/Vbase,'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueV2/Vbase,'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('V2, pu')
    plt.title('True Q-V Curve and Linearization Curve')
    plt.legend()
    plt.savefig('True_QV_Curve_and_Linearization_Curve.png')
    
    plt.figure(4)
    plt.plot(Q12/Sbase, lznDel2,'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueDel2,'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('Delta2, degrees')
    plt.title('True Q-Del Curve and Linearization Curve')
    plt.legend()
    plt.savefig('True_Q-Del_Curve_and_Linearization_Curve.png')
    
    
    solns['trueV2'] = trueV2
    solns['trueDel2'] = trueDel2 
    solns['lznV2'] = lznV2 
    solns['lznDel2'] = lznDel2
    
    return Q12, solns # end of make QV curve


def makePVcurve_3ph(sweep_lb, sweep_ub, Sbase, Vbase, R12, X12, B12, V1):
    # all inputs are NOT in per unit, sweeps are scalars
    numPts = 20
    P12 = np.linspace(sweep_lb, sweep_ub, numPts)
    Q12 = P12*m.tan(m.acos(.9)) # Q=P*tan(acos(0.9)), note P12 and Q12 are vectors of pow vals over time
    R12_diag = np.array([(R12[0][0]), (R12[1][1]), (R12[2][2])])
    X12_diag = np.array([(X12[0][0]), (X12[1][1]), (X12[2][2])])
    #V1_trans = np.transpose(V1)
    V1_trans = (Vbase * np.ones((1,3)))[0]
    trueV2 = np.zeros((3, numPts))
    trueDel2 = np.zeros((3, numPts))
    lznV2 = np.zeros((3, numPts))
    lznDel2 = np.zeros((3, numPts))
    solns = {}
    
    print('HI, R12=',R12)
    print('HI, X12=',X12)
    print('HI, V1=',V1) 
    print('HI, P12/Sbase=',P12[3]/Sbase)
    print('HI, Q12/Sbase-',Q12[3]/Sbase) 

    for i in range(len(P12)): 
        # calculating true curve
        a, b = solveFwdBwdSweep_2bus_3ph(R12, X12, B12, V1, P12[i], Q12[i]) # all in not-pu
        # a is V2, b is del2, V2 is NOT in pu
        trueV2[0][i], trueV2[1][i], trueV2[2][i] = a[0][0], a[1][0], a[2][0]
        trueDel2[0][i], trueDel2[1][i], trueDel2[2][i] = b[0][0], b[1][0], b[2][0]
        
        # calculating linearization curve
        V2sq = (V1_trans**2) - (2*P12[i]*R12_diag) - (2*Q12[i]*X12_diag)
        V2 = V2sq**(1/2)
        lznV2[0][i], lznV2[1][i], lznV2[2][i] = V2[0], V2[1], V2[2]
        delta2 = np.array([0,-120*(m.pi/180),120*(m.pi/180)]) - (((P12[i]*X12_diag)-(Q12[i]*R12_diag))/(V1_trans*V2))
        delta_deg = (180/m.pi)*delta2
        lznDel2[0][i], lznDel2[1][i], lznDel2[2][i] = delta_deg[0], delta_deg[1], delta_deg[2]
    
    plt.figure(figsize = (30,30))
    plt.subplot(431)
    plt.plot(P12/Sbase, lznV2[0]/Vbase,'r', label = 'linearization')
    plt.plot(P12/Sbase, trueV2[0]/Vbase,'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('V2, pu')
    plt.title('True P-V and Linearization Curves: Phase A')
    plt.legend()
    plt.subplot(432)
    plt.plot(P12/Sbase, lznV2[1]/Vbase,'r', label = 'linearization')
    plt.plot(P12/Sbase, trueV2[1]/Vbase,'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('V2, pu')
    plt.title('True P-V and Linearization Curves: Phase B')
    plt.legend()
    plt.subplot(433)
    plt.plot(P12/Sbase, lznV2[2]/Vbase,'r', label = 'linearization')
    plt.plot(P12/Sbase, trueV2[2]/Vbase,'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('V2, pu')
    plt.title('True P-V and Linearization Curves: Phase C')
    plt.legend()
    
    plt.subplot(434)
    plt.plot(P12/Sbase, lznDel2[0],'r', label = 'linearization')
    plt.plot(P12/Sbase, trueDel2[0],'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('Delta2, degrees')
    plt.title('True P-Del and Linearization Curves: Phase A')
    plt.legend()
    plt.subplot(435)
    plt.plot(P12/Sbase, lznDel2[1],'r', label = 'linearization')
    plt.plot(P12/Sbase, trueDel2[1],'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('Delta2, degrees')
    plt.title('True P-Del and Linearization Curves: Phase B')
    plt.legend()
    plt.subplot(436)
    plt.plot(P12/Sbase, lznDel2[2],'r', label = 'linearization')
    plt.plot(P12/Sbase, trueDel2[2],'b', label = 'true')
    plt.xlabel('P12, kW')
    plt.ylabel('Delta2, degrees')
    plt.title('True P-Del and Linearization Curves: Phase C')
    plt.legend()
   
    plt.savefig('True_PV_P-Del_Curve_and_Linearization_Curve.png')
    
    solns['trueV2'] = trueV2
    solns['trueDel2'] = trueDel2 
    solns['lznV2'] = lznV2 
    solns['lznDel2'] = lznDel2
    
    return P12, solns # end of makePVcurve
    
    
def makeQVcurve_3ph(Sweep_lb, Sweep_ub, Sbase, Vbase, R12, X12, B12, V1):
    # All units not in pu, sweep_lb/ub are scalars, R and X are 3x3 matrices
    numPts = 20
    Q12 = np.linspace(Sweep_lb, Sweep_ub, numPts)
    P12 = Q12/m.tan(m.acos(.9)) # P=Q/tan(acos(0.9)), P12 is a vector the size of Q12
    R12_diag = np.array([(R12[0][0]), (R12[1][1]), (R12[2][2])])
    X12_diag = np.array([(X12[0][0]), (X12[1][1]), (X12[2][2])])
    #V1_trans = np.transpose(V1)
    V1_trans = (Vbase * np.ones((1,3)))[0]
    trueV2 = np.zeros((3, numPts))
    trueDel2 = np.zeros((3, numPts))
    lznV2 = np.zeros((3, numPts))
    lznDel2 = np.zeros((3, numPts))
    solns = {}
    
    for i in range(len(Q12)):
        # calculating true curve
        a, b = solveFwdBwdSweep_2bus_3ph(R12, X12, B12, V1, P12[i], Q12[i])  # all in not-pu
        # a is V2, b is del2, V2 is NOT in pu
        trueV2[0][i], trueV2[1][i], trueV2[2][i] = a[0][0], a[1][0], a[2][0]
        trueDel2[0][i], trueDel2[1][i], trueDel2[2][i] = b[0][0], b[1][0], b[2][0]
        
        # calculating linearization curve
        V2sq = (V1_trans**2) - (2*P12[i]*R12_diag) - (2*Q12[i]*X12_diag)
        V2 = V2sq**(1/2)
        delta2 = np.array([0,-120*(m.pi/180),120*(m.pi/180)]) - (((P12[i]*X12_diag)-(Q12[i]*R12_diag))/(V1_trans*V2))
        lznV2[0][i], lznV2[1][i], lznV2[2][i] = V2[0], V2[1], V2[2]
        delta_deg = (180/m.pi)*delta2
        lznDel2[0][i], lznDel2[1][i], lznDel2[2][i] = delta_deg[0], delta_deg[1], delta_deg[2]
    
        
    plt.figure(figsize = (30,30))
    plt.subplot(437)
    plt.plot(Q12/Sbase, lznV2[0]/Vbase,'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueV2[0]/Vbase,'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('V2, pu')
    plt.title('True Q-V and Linearization Curves: Phase A')
    plt.legend()
    plt.subplot(438)
    plt.plot(Q12/Sbase, lznV2[1]/Vbase,'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueV2[1]/Vbase,'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('V2, pu')
    plt.title('True Q-V and Linearization Curves: Phase B')
    plt.legend()
    plt.subplot(439)
    plt.plot(Q12/Sbase, lznV2[2]/Vbase,'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueV2[2]/Vbase,'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('V2, pu')
    plt.title('True Q-V and Linearization Curves: Phase C')
    plt.legend()
    
    plt.subplot(4,3,10)
    plt.plot(Q12/Sbase, lznDel2[0],'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueDel2[0],'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('Delta2, degrees')
    plt.title('True Q-Del and Linearization Curves: Phase A')
    plt.legend()
    plt.subplot(4,3,11)
    plt.plot(Q12/Sbase, lznDel2[1],'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueDel2[1],'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('Delta2, degrees')
    plt.title('True Q-Del and Linearization Curves: Phase B')
    plt.legend()
    plt.subplot(4,3,12)
    plt.plot(Q12/Sbase, lznDel2[2],'r', label = 'linearization')
    plt.plot(Q12/Sbase, trueDel2[2],'b', label = 'true')
    plt.xlabel('Q12, kVAR')
    plt.ylabel('Delta2, degrees')
    plt.title('True Q-Del and Linearization Curves: Phase C')
    plt.legend()
    
    plt.savefig('True_QV_Q-Del_Curve_and_Linearization_Curve.png')
    
    solns['trueV2'] = trueV2
    solns['trueDel2'] = trueDel2 
    solns['lznV2'] = lznV2 
    solns['lznDel2'] = lznDel2
    
    return Q12, solns # end of make QV curve


def computeLznItvl(x, fx_lzn, fx_true):
    # compute lzn itvl:
    err = abs(fx_lzn - fx_true)
    err_max = max(err)
    # compute slopes everywhere on true curve:
    r = [(fx_true[1] - fx_true[0]) / (x[1] - x[0])]
    
    for i in range(len(fx_true) - 2):
        r += [(fx_true[i + 2] - fx_true[i]) / (x[i + 2] - x[i])]
    
    my_index = len(fx_true) - 2
    r += [(fx_true[my_index] - fx_true[my_index + 1]) / (x[my_index] - x[my_index + 1])]
    slopes_true = r
    slopeInfo = [min(slopes_true), max(slopes_true), np.mean(slopes_true)]
    
    return err_max, slopeInfo


def detLznRange(feeder, Vbase_ll, Sbase, z12, B12, act_locs, load_data, headerpath, substation_name, modelpath, depths):
    Vbase = Vbase_ll/(3**(1 / 2)) 
    Zbase = (Vbase**2)/(Sbase*1000)
    V1a=Vbase*np.cos(0*np.pi/180)+1j*Vbase*np.sin(0*np.pi/180); # angle=0 deg
    V1b=Vbase*np.cos(-120*np.pi/180)+1j*Vbase*np.sin(-120*np.pi/180); # angle=-120 deg
    V1c=Vbase*np.cos(+120*np.pi/180)+1j*Vbase*np.sin(+120*np.pi/180); # angle=+120 deg
    V1=np.transpose([[V1a, V1b, V1c]]); # not pu
    R12 = z12.real # not pu
    X12 = z12.imag
    
    P_lb_results, P_ub_results, Q_lb_results, Q_ub_results = computePQsweep_timesteps(feeder, load_data)
    PQ_bounds = computePQsweep_losses(feeder, act_locs, Sbase, Zbase, P_lb_results, P_ub_results, Q_lb_results, Q_ub_results, headerpath, substation_name, modelpath, depths)
    
    Psweep_lb = PQ_bounds[0] # in pu
    Psweep_ub = PQ_bounds[1] # each is scalar
    Qsweep_lb = PQ_bounds[2]
    Qsweep_ub = PQ_bounds[3]
    pvals, solns1 = makePVcurve_3ph(Psweep_lb*Sbase, Psweep_ub*Sbase, Sbase, Vbase, R12, X12, B12, V1) # all in not-pu
    qvals, solns2 = makeQVcurve_3ph(Qsweep_lb*Sbase, Qsweep_ub*Sbase, Sbase, Vbase, R12, X12, B12, V1)
    
    errVmax1, slope_vp = computeLznItvl(pvals / Sbase, solns1['lznV2'][0] / Vbase, solns1['trueV2'][0] / Vbase)
    errDelmax1,slope_delp = computeLznItvl(pvals / Sbase, solns1['lznDel2'][0], solns1['trueDel2'][0])
    
    errVmax2, slope_vq = computeLznItvl(qvals / Sbase, solns2['lznV2'][0] / Vbase, solns2['trueV2'][0] / Vbase)
    errDelmax2, slope_delq = computeLznItvl(qvals / Sbase, solns2['lznDel2'][0], solns2['trueDel2'][0])
    
    true_lzn_error_max = [errVmax1, errVmax2, errDelmax1, errDelmax2]
    slope_across_PQ_sweep = [slope_vp, slope_delp, slope_vq, slope_delq]
    
    return true_lzn_error_max, slope_across_PQ_sweep 