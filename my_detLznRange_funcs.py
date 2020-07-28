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
import my_heatmapSetup_funcs as hm



def compute_line_losses_multiphase(feeder, P_vals, Q_vals, act_locs, lb_mode = True):
    graph = feeder.network
    n = len(graph.nodes) 
    percent_V_drop = .05
    P = P_vals
    Q = Q_vals
    load_data_table = pd.read_csv(loadpath, header = None) # file used to extract table headers (node names) ==> not used to calculate actual loads
    num_cols = len(load_data_table.columns)
    branch_list = segment_network(feeder, substation_name)
    status = 'unsolved'
    
    node_index_map = createNodeIndexMap(feeder)
    headers = []  
    
    # retrieve headers from load data file
    for col_index in range(num_cols)[1:]:
        col = load_data_table[load_data_table.columns[col_index]].tolist()
        my_val = col[0]
        headers += [my_val]
    
    S = np.zeros((3, n))
    S = S.astype('complex128')
    
    # parse through load data headers to find node names
    for i in range(len(P)):
        label = headers[i]
        phase = int(label[len(label):len(label)-2:-1]) - 1
        node_name = label[3:]
        node_name = node_name[len(node_name)-4::-1]
        node_name = node_name[len(node_name)::-1]
        node_name = 'bus_' + node_name
        index = node_index_map[node_name]
        S[phase][index] = complex(P[i], Q[i])
    
    Sact = np.zeros(n)
    
    # set actuators to extreme values for lower or upper bound
    if lb_mode:
        act_cap = -.3
    
    else:
        act_cap = .3
    
    for act in act_locs:
        index = node_index_map[act]
        Sact[index] = act_cap
        
    # Find edge nodes:
    edge_nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            edge_nodes += [node]
    
    # Find impedances from edge nodes to substation
    z_edges_to_sub = np.zeros((3, len(edge_nodes)))
    i = 0
    
    for edge in edge_nodes:
        cur_z = imp.get_total_impedance_from_substation(feeder, edge)
        z_edges_to_sub[0][i] = cur_z['Phase 1'] 
        z_edges_to_sub[1][i] = cur_z['Phase 2'] 
        z_edges_to_sub[2][i] = cur_z['Phase 3']
        i += 1
    
    zmax = max([max(z_edges_to_sub[0]), max(z_edges_to_sub[1]), max(z_edges_to_sub[2])])
    Vest = np.zeros((3, n)) # Vest is proportional to distance from substaion.
    
    while percent_V_drop >= -.05 and status == 'unsolved':
        loop_status = 'unbroken'
        
        # estimate voltages at edge nodes
        for i in range(len(edge_nodes)):
            cur_edge = edge_nodes[i]
            cur_index = node_index_map[cur_edge]
            cur_za, cur_zb , cur_zc = z_edges_to_sub[0][i], z_edges_to_sub[1][i], z_edges_to_sub[2][i]
            va = 1 - (percent_V_drop * (cur_za / zmax)) # per unit
            vb = 1 - (percent_V_drop * (cur_zb / zmax))
            vc = 1 - (percent_V_drop * (cur_zc / zmax))
            Vest[0][cur_index] = va 
            Vest[1][cur_index] = vb
            Vest[2][cur_index] = vc 

        V = Vest
        Sloss = np.zeros((n,n)) # will populate 
        Sloss = Sloss.astype('complex128')
        Vforks = {}
        Iforks = {}
        active_branches = []
        
        for branch in branch_list:
            
            for edge in edge_nodes:
                
                if edge in branch:
                    active_branches += [branch]
                    break
                    
            if len(branch) >= 2:       
                edge = branch[len(branch) - 1: len(branch) - 2: -1]
            
            else:
                edge = [branch[0]]
                
            Iforks[edge[0]] = []
            Vforks[edge[0]] = []

        for _ in range(len(branch_list)):
            if loop_status == 'broken':
                break
                
            cur_branch = active_branches[0]
            print(cur_branch)
            
            if len(cur_branch) >= 2:
                edge = cur_branch[len(cur_branch) - 1: len(cur_branch) - 2: -1][0]
                
            else:
                edge = cur_branch[0]
                
            pred_node = list(graph.predecessors(edge))
            cur_index = node_index_map[edge]
            pred_index = node_index_map[pred_node[0]]
            impedance = feeder.network.get_edge_data(pred_node[0], edge, default=None)['connector']
            z_3by3 = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3))

            ia_conj = (S[0][cur_index] + Sact[cur_index]) / V[0][cur_index]
            ib_conj = (S[1][cur_index] + Sact[cur_index]) / V[1][cur_index]
            ic_conj = (S[2][cur_index] + Sact[cur_index]) / V[2][cur_index]
            # take complex conjugate:
            ia_load = complex(ia_conj.real, (-1) * ia_conj.imag)
            ib_load = complex(ib_conj.real, (-1) * ib_conj.imag)
            ic_load = complex(ic_conj.real, (-1) * ic_conj.imag)
            
            ia = ia_load + sum([i[0] for i in Iforks[edge]])
            ib = ib_load + sum([i[1] for i in Iforks[edge]])
            ic = ic_load + sum([i[2] for i in Iforks[edge]])
            
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

            while len(cur_node_children) == 1 and cur_node != substation_name:
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
                z_3by3 = impedance.Z if isinstance(impedance, setup_nx.line) else np.zeros((3,3))
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

            if cur_node == substation_name:
                Seq = sum(sum(Sloss)) + sum(sum(S)) + sum(Sact)
                Peq = Seq.real
                Qeq = Seq.imag
                return Peq, Qeq
            
            if len(cur_node_children) > 1:
                V_drop = np.dot(z_3by3, np.array([ia_prev, ib_prev, ic_prev]))
                Vfork = [V[0][child_index] + V_drop[0]]
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
                            
                            if loop_status == 'broken':
                                break
                                
                            for v_inner in val:
                                checkVclose = [((v_outer[0] - v_inner[0]) / v_outer[0]) < 1] # check not more than 10% different
                                checkVclose += [((v_outer[1] - v_inner[1]) / v_outer[1]) < 1]
                                checkVclose += [((v_outer[2] - v_inner[2]) / v_outer[2]) < 1]

                                if False in checkVclose:
                                    percent_V_drop -= .001
                                    loop_status = 'broken'
                                    break

                        if loop_status == 'unbroken':
                            key_index = node_index_map[key]
                            V[0][key_index] = sum([v[0] for v in val]) / len(val)
                            V[1][key_index] = sum([v[1] for v in val]) / len(val)
                            V[2][key_index] = sum([v[2] for v in val]) / len(val)
                            
                            for branch in branch_list:

                                if key in branch:
                                    active_branches += [branch]
                                    Vforks[key] = []
                                    print(key)
                                    break

    if percent_V_drop < -.05:
        print('Voltages did not converge at forks')
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
    data_table = pd.read_csv(load_data, header = None)
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
    
    Psweep_lb = min(run_sum_P)
    Psweep_ub = max(run_sum_P)
    Qsweep_lb = min(run_sum_Q)
    Qsweep_ub = max(run_sum_Q)
    
    time_P_lb, P_lb_rowP, P_lb_rowQ = find_timestep(Psweep_lb, run_sum_P, Qstart_index, data_table) 
    time_P_ub, P_ub_rowP, P_ub_rowQ = find_timestep(Psweep_ub, run_sum_P, Qstart_index, data_table)
    time_Q_lb, Q_lb_rowP, Q_lb_rowQ = find_timestep(Qsweep_lb, run_sum_Q, Qstart_index, data_table)
    time_Q_ub, Q_ub_rowP, Q_ub_rowQ = find_timestep(Qsweep_ub, run_sum_Q, Qstart_index, data_table)
    
    P_lb_results = [time_P_lb, P_lb_rowP, P_lb_rowQ]
    P_ub_results = [time_P_ub, P_ub_rowP, P_ub_rowQ]
    Q_lb_results = [time_Q_lb, Q_lb_rowP, Q_lb_rowQ]
    Q_ub_results = [time_Q_ub, Q_ub_rowP, Q_ub_rowQ]
    
    print('Timesteps for Extreme Load Values:')
    print('P_lb = ',time_P_lb)
    print('P_ub = ',time_P_ub)
    print('Q_lb = ',time_Q_lb)
    print('Q_ub = ',time_Q_ub)
    
    return P_lb_results, P_ub_results, Q_lb_results, Q_ub_results # end of computing line losses


def computePQsweep_losses(feeder, act_locs, P_lb_results, P_ub_results, Q_lb_results, Q_ub_results):
    time_P_lb, P_lb_rowP, P_lb_rowQ = P_lb_results[0], P_lb_results[1], P_lb_results[2]
    time_P_ub, P_ub_rowP, P_ub_rowQ = P_ub_results[0], P_ub_results[1], P_ub_results[2]
    time_Q_lb, Q_lb_rowP, Q_lb_rowQ = Q_lb_results[0], Q_lb_results[1], Q_lb_results[2]
    time_Q_ub, Q_ub_rowP, Q_ub_rowQ = Q_ub_results[0], Q_ub_results[1], Q_ub_results[2]
    
    # accounting for line losses:
    Psweep_lb, _ = compute_line_losses_multiphase(feeder, P_lb_rowP, P_lb_rowQ, act_locs, lb_mode = True)
    Psweep_ub, _ = compute_line_losses_multiphase(feeder, P_ub_rowP, P_ub_rowQ, act_locs, lb_mode = False)
    _, Qsweep_lb = compute_line_losses_multiphase(feeder, Q_lb_rowP, Q_lb_rowQ, act_locs, lb_mode = True)
    _, Qsweep_ub = compute_line_losses_multiphase(feeder, Q_ub_rowP, Q_ub_rowQ, act_locs, lb_mode = False)
    
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
def solveFwdBwdSweep_2bus_3ph(R12, X12, Vs, P2, Q2):
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
    I12 = np.conj(np.divide(np.transpose(S2),V2[k]))
    #print('I12=',I12)
    #print(np.concatenate((Vconv[2*k], Vconv[2*k+1])))
    
    '''Iterative Part'''
    
    while any(node >= tol for node in np.concatenate((Vconv[2*k], Vconv[2*k+1]))):  # break when all nodes less than tol
        k += 1  # new iteration
         # Fwd sweep
        V1=np.append(V1,[V1[k-1]],axis=0)# same as prev iter ZERO?
        V2=np.append(V2,Vs - np.dot(I12,z12),axis=0)
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
        I12 = np.conj(np.divide(np.transpose(S2),V2[k]))

        if len(Vconv) > 30:
            print('Didnt converge')
        break  # break out of loop of iter too many times

    '''Output Results'''
    #print('~~~~~~~ PF Results: ')
    V1soln=V1[-1]
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
        a, b = solveFwdBwdSweep_2bus(R12, X12, V1, P12[i], Q12)
        trueV2[i] = a
        trueDel2[i] = b
        V2sq = (V1**2) - (2*R12*P12[i]) - (2*X12*Q12)
        V2 = V2sq**(1/2)
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


def makePVcurve_3ph(sweep_lb, sweep_ub, Sbase, Vbase, R12, X12, V1):
    numPts = 20
    P12 = Sbase * np.linspace(sweep_lb, sweep_ub, numPts)
    Q12pu = m.tan(m.acos(.9))
    Q12 = Q12pu * Sbase
    Zbase = (Vbase**2)/Sbase
    R12_diag = Zbase * np.array([[(R12[0][0]), (R12[1][1]), (R12[2][2])]])
    X12_diag = Zbase * np.array([[(X12[0][0]), (X12[1][1]), (X12[2][2])]])
    V1_trans = Vbase * np.transpose(V1)
    trueV2 = np.zeros((3, numPts))
    trueDel2 = np.zeros((3, numPts))
    lznV2 = np.zeros((3, numPts))
    lznDel2 = np.zeros((3, numPts))
    solns = {}
    
    for i in range(len(P12)):
        a, b = solveFwdBwdSweep_2bus_3ph(R12, X12, V1, P12[i]/Sbase, Q12/Sbase)
        a = a*Vbase
        trueV2[0][i], trueV2[1][i], trueV2[2][i] = a[0][0], a[1][0], a[2][0]
        trueDel2[0][i], trueDel2[1][i], trueDel2[2][i] = b[0][0], b[1][0], b[2][0]
        V2sq = (V1_trans**2) - (2*P12[i]*R12_diag) - (2*Q12*X12_diag)
        V2 = V2sq**(1/2)
        delta2 = -1 * (((P12[i]*X12_diag)-(Q12*R12_diag))/(V1_trans*V2))
        lznV2[0][i], lznV2[1][i], lznV2[2][i] = V2[0][0], V2[0][1], V2[0][2]
        delta_deg = (180/m.pi)*delta2
        lznDel2[0][i], lznDel2[1][i], lznDel2[2][i] = delta_deg[0][0], delta_deg[0][1], delta_deg[0][2]
    
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
    
    
def makeQVcurve_3ph(Sweep_lb, Sweep_ub, Sbase, Vbase, R12, X12, V1):
    numPts = 20
    Q12 = Sbase * np.linspace(Sweep_lb, Sweep_ub, numPts)
    P12pu = m.tan(m.acos(0.9))  
    P12 = P12pu * Sbase
    Zbase = (Vbase**2)/Sbase
    R12_diag = Zbase * np.array([[(R12[0][0]), (R12[1][1]), (R12[2][2])]])
    X12_diag = Zbase * np.array([[(X12[0][0]), (X12[1][1]), (X12[2][2])]])
    V1_trans = Vbase * np.transpose(V1)
    trueV2 = np.zeros((3, numPts))
    trueDel2 = np.zeros((3, numPts))
    lznV2 = np.zeros((3, numPts))
    lznDel2 = np.zeros((3, numPts))
    solns = {}
    
    for i in range(len(Q12)):
        a, b = solveFwdBwdSweep_2bus_3ph(R12, X12, V1, P12/Sbase, Q12[i]/Sbase)
        a = a*Vbase
        trueV2[0][i], trueV2[1][i], trueV2[2][i] = a[0][0], a[1][0], a[2][0]
        trueDel2[0][i], trueDel2[1][i], trueDel2[2][i] = b[0][0], b[1][0], b[2][0]
        V2sq = (V1_trans**2) - (2*P12*R12_diag) - (2*Q12[i]*X12_diag)
        V2 = V2sq**(1/2)
        delta2 = -1 * (((P12*X12_diag)-(Q12[i]*R12_diag))/(V1_trans*V2))
        lznV2[0][i], lznV2[1][i], lznV2[2][i] = V2[0][0], V2[0][1], V2[0][2]
        print(delta2)
        delta_deg = (180/m.pi)*delta2
        lznDel2[0][i], lznDel2[1][i], lznDel2[2][i] = delta_deg[0][0], delta_deg[0][1], delta_deg[0][2]
    
        
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


def detLznRange(feeder, Vbase_ll, Sbase, z12, act_locs):
    Vbase = Vbase_ll/(3**(1 / 2)) # not pu
    V1 = Vbase * np.ones((3,1)) # slack bus, not pu
    R12 = z12.real
    X12 = z12.imag
    
    PQ_bounds, PQ_extremes_time_steps = computePQsweep(feeder, load_data, act_locs)
    Psweep_lb = PQ_bounds[0]
    Psweep_ub = PQ_bounds[1]
    Qsweep_lb = PQ_bounds[2]
    Qsweep_ub = PQ_bounds[3]
    
    pvals, solns1 = makePVcurve_3ph(Psweep_lb, Psweep_ub, Sbase, Vbase, R12, X12, V1)
    qvals, solns2 = makeQVcurve_3ph(Qsweep_lb, Qsweep_ub, Sbase, Vbase, R12, X12, V1)
    
    errVmax1, slope_vp = computeLznItvl(pvals / Sbase, solns1['lznV2'] / Vbase, solns1['trueV2'] / Vbase)
    errDelmax1,slope_delp = computeLznItvl(pvals / Sbase, solns1['lznDel2'], solns1['trueDel2'])
    
    errVmax2, slope_vq = computeLznItvl(qvals / Sbase, solns2['lznV2'] / Vbase, solns2['trueV2'] / Vbase)
    errDelmax2, slope_delq = computeLznItvl(qvals / Sbase, solns2['lznDel2'], solns2['trueDel2'])
    
    true_lzn_error_max = [errVmax1, errVmax2, errDelmax1, errDelmax2]
    slope_across_PQ_sweep = [slope_vp, slope_delp, slope_vq, slope_delq]
    
    return true_lzn_error_max, slope_across_PQ_sweep 