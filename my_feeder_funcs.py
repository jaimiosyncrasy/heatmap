# Hi
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

print(type(pd))

def feeder_init(modelpath,loadfolder,loadpath,timesteps,Vbase_ll,Sbase,depths,leaves):

    modeldata = pd.ExcelFile(modelpath)
    actpath = loadpath
    
    # set dummy values for undefined variables
    date = datetime.datetime.now()
    month = date.month
    day = date.day
    hour = date.hour
    minute = date.minute
    #timestepcur = hour*60+minute
    timestepcur = 0
    
    Psat_nodes = []
    Qsat_nodes = []
    
    refphasor = np.ones((3,2))
    refphasor[:,0]=1
    refphasor[:,1]=[0,4*np.pi/3,2*np.pi/3]
    
    PVforecast = 0
    
    #get feeder
    fin_feeder = feeder(modelpath,loadfolder,loadpath,actpath,timesteps,timestepcur,
                         Vbase_ll,Sbase,refphasor,Psat_nodes,Qsat_nodes,PVforecast, depths, leaves)    
    return fin_feeder

def nx_plot(feeder):
    #plt.title('13_node_test')
    #pos = nx.nx_agraph.graphviz_layout(feeder.network, prog='dot')
    #pos = nx.spectral_layout(feeder.network)
    #nx.draw(feeder.network, pos, with_labels=True, arrows=True)
    #plt.savefig('13_node_test.png')
    
    ### THESE ARE ADDITIONS ONLY FOR THE PL0001 FEEDER
    #feeder.network.add_edge('bus_N_300062298', 'bus_N_L_22666_sec')
    #feeder.network.add_edge('bus_N_300062062', 'bus_N_L_21316_sec')
    #feeder.network.add_edge('bus_N_300062332', 'bus_N_L_22077_sec')
    #feeder.network.add_edge('bus_N_300053280', 'bus_N_L_52586_sec')
    #feeder.network.add_edge('bus_N_300006691', 'bus_N_L_38426_sec')
    
    ### THESE ARE ADDITIONS ONLY FOR THE AL0001 FEEDER
    #feeder.network.add_edge('bus_N_900081044', 'bus_N_L_87632_sec')
    #feeder.network.add_edge('bus_N_900059208', 'bus_N_L_46793_sec')
    #feeder.network.add_edge('bus_N_900060818', 'bus_N_L_17532_sec')
    #feeder.network.add_edge('bus_N_900076500', 'bus_N_L_108238_sec')
    #feeder.network.add_edge('bus_N_900059040', 'bus_N_L_111610_sec')
    #feeder.network.add_edge('bus_N_900059556', 'bus_N_L_147411_sec')
    #feeder.network.add_edge('bus_N_900047073', 'bus_N_L_6860_sec')
    #feeder.network.add_edge('bus_N_900054239', 'bus_N_L_9709_sec')
    #feeder.network.add_edge('bus_N_900056678', 'bus_N_L_110483_sec')
    #feeder.network.add_edge('bus_N_900059242', 'bus_N_L_132901_sec')
    #feeder.network.add_edge('bus_N_900019520', 'bus_N_L_40688_sec')
    #feeder.network.add_edge('bus_N_900056622', 'bus_N_L_138443_sec')
    #feeder.network.add_edge('bus_N_900059203', 'bus_N_L_108120_sec')
    #feeder.network.add_edge('bus_N_900080808', 'bus_N_L_116563_sec')
    #feeder.network.add_edge('bus_N_900008961', 'bus_N_L_113827_sec')
    
    #nx.nx_pydot.write_dot(feeder.network, '13_node_test.dot')
    return

def make_graph(feeder, file_name):
    #file_name = string specifying name of dot file created when make_graph() is called
    #feeder = initialized feeder object from feeder_init()
    #for node in feeder.network.nodes:
    graph = feeder.network
    nodes = graph.nodes
    edges = graph.edges
    
    for n in nodes:
        graph.nodes[n]['fontsize'] = 15
        graph.nodes[n]['label'] = ''
        graph.nodes[n]['shape'] = 'circle'
        graph.nodes[n]['width'] = .25
        graph.nodes[n]['xlabel'] = n[4:]
        
    for e in edges:
        graph.edges[e]['arrowhead'] = 'none'
        graph.edges[e]['color'] = 'grey'
    
    nx.nx_pydot.write_dot(graph, file_name)
    render('dot', 'png', file_name)  
    #nx.draw_planar(graph,with_labels = True, alpha=0.8) #NEW FUNCTION
    #nx.draw_circular(graph,with_labels = True, alpha=0.8) #NEW FUNCTION
    #nx.draw_shell(graph,with_labels = True, alpha=0.8) #NEW FUNCTION
    #nx.draw_kamada_kawai(graph,with_labels = True, alpha=0.8) #NEW FUNCTION
    #nx.draw_spectral(graph,with_labels = True, alpha=0.8) #NEW FUNCTION
    #nx.draw_spring(graph,with_labels = True, alpha=0.8) #NEW FUNCTION
    return

def clear_graph(feeder):
    # returns feeder graph to original state (removes any node colors/changes to shape or size)
    graph = feeder.network
    for n in graph.nodes:
        graph.nodes[n]['style'] = 'filled'
        graph.nodes[n]['fillcolor'] = 'white'
        graph.nodes[n]['color'] = 'black'
        graph.nodes[n]['fontcolor'] = 'black'
        graph.nodes[n]['fontsize'] = 15
        graph.nodes[n]['label'] = ''
        graph.nodes[n]['shape'] = 'circle'
        graph.nodes[n]['width'] = .25
    return

# see reference for layouts: https://networkx.github.io/documentation/stable/reference/drawing.html
# see img: https://i.stack.imgur.com/3ZRVT.png