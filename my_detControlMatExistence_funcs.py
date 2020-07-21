import importlib
import setup_nx # your own module, setup.nx.py
import numpy as np
import scipy.linalg as LA
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
import my_detLznRange_funcs as lzn
import my_heatmapSetup_funcs as hm


def assignF(Fp,Fq,indicMat):
    n=int(len(indicMat)/6) # indicMat is 6n x 6n
    
    Hblock1=np.zeros((3*n,3*n))   
    ridx,colidx=np.nonzero(indicMat[0:3*n,0:3*n]) # python indexing goes first to (last-1) 
    Hblock1[ridx,colidx]=Fq

    Hblock2=np.zeros((3*n,3*n))   
    ridx,colidx=np.nonzero(indicMat[0:3*n,0:3*n]) 
    Hblock2[ridx,colidx]=Fp
  
    upper=np.concatenate((Hblock1, np.zeros((3*n,3*n))),axis=1)
    lower=np.concatenate((np.zeros((3*n,3*n)),Hblock2),axis=1)
    F=np.concatenate((upper,lower))
    
    print(F)
    #print("Size of F=",F.shape)
    return F

# version 1
def computeFParamSpace_v1(feeder, act_locs, perf_nodes):
    # Compute (Fp,Fq) ranges as a func of impedance paths between act nodes, perf nodes, and substation
    gamma=[0.08, 0.35]
    Fplist,Fqlist= (np.empty((0,1)) for i in range(2))
    for i in len(act_locs):
        Fqlist[i]=i
        Fplist[i]=i
        
    Fq_lb=np.min(Fqlist)
    Fp_lb=np.min(Fplist)
    return Fq_lb,Fp_lb

#version 2
def computeFParamSpace_v2(feeder, act_locs, perf_nodes,R,X,depths,node_index_Map):
    # Compute (Fp,Fq) ranges as a func of impedance paths between act nodes, perf nodes, and substation
    gamma=[0.08, 0.35] # scaling, according to units of Q-V and P-delta loops
    c=[1, 1]
    sensEst_dvdq,sensEst_deldp= (np.empty((0,1)) for i in range(2))
    
    for act in act_locs: # for each (perf,act) pair
        print(node_index_Map)
        act_idx=node_index_Map[act] # index of act node
        print('act idx=',act_idx)

        perf_idx=act_idx #index of act node
        perf=perf_nodes[0] # temp
        # indexing  a[start:stop] --> items start through stop-1
        print(R)
        print(len(R))
        Zgood=np.empty((3,3))
        Zgood=R[act_idx:(act_idx+3),perf_idx:(perf_idx+3)]+X[act_idx:(act_idx+3),perf_idx:(perf_idx+3)]*1j # 3x3 matrix
        #Zgood=R[act_idx:(act_idx+3),perf_idx:(perf_idx+3)]
        print('Zgood=\n',Zgood)
        
        Z_toSubst=imp.get_total_impedance_from_substation(feeder,act,depths)
        Z_actperf=imp.get_total_impedance_between_two_buses(feeder,act,perf,depths)
        
        # until this func returns 3x3 array...
#         Z_toSubst=np.zeros([3,3],dtype=complex)
#         D=imp.get_total_impedance_from_substation(feeder,act,depths)
#         mylist=list(D.values())
#         print(mylist[0])
#         print(Z_toSubst)
#         Z_toSubst[0,0]=mylist[0]
#         Z_toSubst[1,1]=mylist[1]
#         Z_toSubst[2,2]=mylist[2]
#         np.zeros((5,), dtype=int)

#         Z_actperf=np.zeros([3,3],dtype=complex)
#         D=imp.get_total_impedance_between_two_buses(feeder,act,perf,depths)
#         mylist=list(D.values())
#         Z_actperf[0,0]=mylist[0]
#         Z_actperf[1,1]=mylist[1]
#         Z_actperf[2,2]=mylist[2]
        
        #Zgood and Z_tosubts should have same value, but they dont...
        
        print('Z to subst=',Z_toSubst)
        Zbad1=np.subtract(Z_toSubst,Zgood) # >=0
        Zbad2=np.subtract(Z_actperf,Zbad1) # >=0
        #Zbad2=imp.get_total_impedance_between_two_buses(feeder,act,perf,depths) - Zbad1
        print('Zbad1=',Zbad1)
        print('Zbad2=',Zbad2)
#         der1=1-np.divide(c[1]*Zbad1,Zgood+Zbad1) # derating for actuator not colocated
#         der2=1-np.divide(c[2]*Zbad2,Zgood+Zbad2) # derating for perf node not on samepath to substation as act
#         print(der1)
#         print(der2)
#         print(mainPath)
#         sensEst_dvdq[i]=der1*der2*gamma[1]*mainPath # each element is 3x3 matrix
#         sensEst_deldp[i]=der1*der2*gamma[2]*mainPath

#for my_buses,list1, list2,list3 in zip(my_buses,list1,list2,list3):
 #   print(my_buses,list1, list2,list3)

    #Fq_lb=np.min(sensEst_dvdq)
    #Fp_lb=np.min(sensEst_deldp)
    Fq_lb=1
    Fp_lb=2
    return Fq_lb,Fp_lb

#   Compute (Fp,Fq) ranges as a func 
# Zmag=sqrt(X.^2+R.^2);
# gamma=[0.08 0.35]; % scaling factor, different for Fp vs. Fq
# 
# % For one config, may have several act-perf node pairs
# actLoc_vec=[4] % node 3 and 4
# perfLoc_vec=[3]; % node 2 and 6
# Fp=[]; Fq=[];
#  for i=1:numAct
#      actLoc=actLoc_vec(i);
#      perfLoc=perfLoc_vec(i);
#      mainPath=Zmag(perfLoc,perfLoc); % path from perf node to substation z01+z12
#      % deratePath written below is not generalizable to any config. Rewrite
#      % to be the path between perf node and act node
#      deratePath=Zmag(actloc,actloc)-Zmag(perfLoc,perfLoc); % path between perf node and act node
#      derating=(1-min(deratePath,mainPath)/mainPath) % due to actuator not being co0located with perf node
#      Fp(i)=derating*gamma(1)*mainPath
#      Fq(i)=derating*gamma(2)*mainPath
#  end
#  Fp_ub=min(Fp); % range of Fp and Fq chosen as the smallest among all act-perf pairs, smallest because we set all kgains equal so better to be under than over aggressive
#  Fq_ub=min(Fq); 
#  
#  % Fq=0.6, Fp=0.1


def detControlMatExistence(A, B, indicMat):
#def detControlMatExistence(feeder, act_locs, perf_nodes,A,B,R,X,indicMat):
    n=int(len(indicMat)/6) # indicMat is 6n x 6n

# Compute good (Fp,Fq) sample space
    # Fq_lb,Fp_lb=computeFParamSpace_v1(feeder, act_locs, perf_nodes)
    # NOTE: need modify detControlMatExistence input params to pass stuff into computeFparamSpace
    
    #Fq_range=np.arange(-2, 4, 0.2).tolist()
    #Fp_range=np.arange(-0.1, 0.1, 0.015).tolist()
    numSamp=10
    #Fq_range=np.arange(-Fq_lb, Fq_lb, 2*Fq_lb/numSamp).tolist()
    #Fp_range=np.arange(-Fp_lb,Fp_lb, 2*Fp_lb/numSamp).tolist()
    Fq_range=np.arange(-1, 1, 0.5).tolist()
    Fp_range=np.arange(-2, 2, 0.5).tolist()
    
# Initialize arrays, will be populated in loops
    feas=False # boolean
    numfeas,myCosts= (np.empty((0,1)) for i in range(2))
    feasFs,myFbases=(np.empty((0,2)) for i in range(2)) # 0 for dim to concatenate on, 2 for length of vectors being concatenated
    dataFull, data_zero = (np.array([]) for i in range(2)) # for each config
    CLmat=np.empty((6*n,6*n))

    for Fq in Fq_range:
        for Fp in Fp_range:
            if not(Fq==0 or Fp==0): # skip iteration if either zero
                #print("(Fq,Fp)=",Fq,",",Fp,")")
                F=assignF(Fp,Fq,indicMat)
                CLmat=A-np.dot(B,F) # CLmat=A-BF
                eigs,evecs=LA.eig(CLmat) # closed loop eigenvalues
                eigMags=np.absolute(eigs)
              #  if all(np.absolute(eigs)<=np.ones(6*n,1)):
                if all(np.around(eigMags,decimals=6)<=1): 
                    feas=True
                    feasFs=np.append(feasFs,[[Fp, Fq]],axis=0)
                val=np.sum(eigMags[np.where(eigMags > 1)])
                myCosts=np.append(myCosts,[[val]],axis=0) # temp
                myFbases=np.append(myFbases,[[Fp, Fq]],axis=0)

    if feas==True:
        print("Config feasible!")
        numfeas=np.append(numfeas,[[len(feasFs)]],axis=0) # number of rows
    else:
        feas=False
        print("No F found for config --> non-conclusive")
        numfeas=np.append(numfeas,[[0]],axis=0)

    dataFull=np.concatenate((myFbases,myCosts),axis=1) # [Fp Fq myCost]  
    print("[Fp,Fq,myCost]=\n",dataFull)
    numTried=len(dataFull) # number of rows
    num_act=np.count_nonzero(indicMat)/2
       
    # return feas,feasFs,num_act,numfeas,numTried
    return feas,feasFs,numfeas,numTried,num_act