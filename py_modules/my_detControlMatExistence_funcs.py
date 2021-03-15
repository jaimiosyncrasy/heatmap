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


def assignF(parmObj,Fp,Fq,indicMat): # algo similar to updateStateSpace

    n=int(len(indicMat)/6) # indicMat has 6n rows for all versions
    Hblock1=np.zeros((3*n,3*n))   
    ridx,colidx=np.nonzero(indicMat[0:3*n,0:3*n]) # python indexing goes first to (last-1)          
    for k in range(len(ridx)):
            Hblock1[ridx[k]][colidx[k]] = Fq
    
    Hblock2=np.zeros((3*n,3*n))   

    if parmObj.get_version()==1: # PBC       
        ridx,colidx=np.nonzero(indicMat[0:3*n,0:3*n]) 
        for k in range(len(ridx)):
                Hblock2[ridx[k]][colidx[k]] = Fp

        upper=np.concatenate((Hblock1, np.zeros((3*n,3*n))),axis=1)
        lower=np.concatenate((np.zeros((3*n,3*n)),Hblock2),axis=1)
        F=np.concatenate((upper,lower))  # indicMat is 6n x 6n, [Fq 0 ; 0 Fp]
    
    else: # Droop:        
        ridx,colidx=np.nonzero(indicMat[3*n+1:,0:3*n]) 
        for k in range(len(ridx)):
                Hblock2[ridx[k]][colidx[k]] = Fp
        F=np.concatenate((Hblock1,Hblock2),axis=0) # indicMat is now 6n x 3n, [Fq Fp]'

    #print("Size of F=",F.shape)
    return F

# version 1, workaround
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

#version 2, correct
def computeFParamSpace_v2(parmObj,feeder, act_locs, perf_nodes,R,X,depths,node_index_Map, file_name):
    # Compute (Fq,Fp) ranges as a func of impedance paths between act nodes, perf nodes, and substation
    if file_name=='13NFbalanced':
        c=np.array([0.412,0.857]) # (q,p) tuned for 13NF based on data from feas configs
    elif file_name=='123NF':
        #c=np.array([0.3,0.45])  # setting on 9/2/20
        #c=np.array([0.5,0.7]) # place_max_coloc getting not-great results from this
        if parmObj.get_version()==1:  
            c=np.array([0.3,0.45])  # setting back to this on 2/28/21
        else: # droop case
            c=np.array([0.7,1.1])  # setting this on 3/7/21
    elif file_name=='PL0001': 
        c=np.array([0.32,0.65])
    else:
        print('error: c not assigned because couldnt find load file in folder')

        
    avgSens_dvdq,avgSens_ddeldp= (np.empty((0,1)) for i in range(2))
    i=0
    for act in act_locs: # for each (perf,act) pair in current config you're evaluating
        #print(node_index_Map)
        perf=perf_nodes[i] #index of act node
        i=i+1 # increment to update perf node
        # indexing  a[start:stop] --> items start through stop-1
        act_idx=node_index_Map[act] # index of act node
        perf_idx=node_index_Map[perf] 

        #print('evaluating act at ',act,', perf at ',perf)

        Zgood=np.empty((3,3))
        Zgood=(R[act_idx*3:(act_idx*3+3),perf_idx*3:(perf_idx*3+3)])/2+(X[act_idx*3:(act_idx*3+3),perf_idx*3:(perf_idx*3+3)])/2*1j # 3x3 matrix        
        Z_toSubst=imp.get_total_impedance_from_substation(feeder,act,depths)
        Z_actperf=imp.get_total_impedance_between_two_buses(feeder,act,perf,depths)

        Zbad1=np.subtract(Z_toSubst,Zgood) # >=0
        Zbad2=np.subtract(Z_actperf,Zbad1) # >=0
        #print('Zbad1=',np.around(Zbad1,2))
        #print('Zbad2=',np.around(Zbad2,2))
        
        # deratings should range 0 to 1 for each 3x3 elements
        der1=np.ones(3)-np.divide(Zbad1,Zgood+Zbad1) # derating for actuator not colocated
        der2=np.ones(3)-np.divide(Zbad2,Zgood+Zbad2) # derating for perf node not on samepath to substation as act
        #print('der1=',der1)
        #print('der2=',der2)
        mainPath=Zgood
        sensEst_dvdq=np.multiply(np.multiply(der1,der2),mainPath)
        sensEst_ddeldp=np.multiply(np.multiply(der1,der2),mainPath)

        avg_dvdq=sensEst_dvdq.mean() # converts 3x3 to scalar complex for each act-perf pair
        avg_ddeldp=sensEst_ddeldp.mean()
        
        if np.isnan(avg_dvdq) or np.isnan(avg_ddeldp):
            print('act=',act)  
            sys.exit('Error: avg_dvdq or avg_ddeldp are NaN')  
            
        avgSens_dvdq=np.append(avgSens_dvdq,avg_dvdq)
        avgSens_ddeldp=np.append(avgSens_ddeldp,avg_ddeldp)
        #print('avgdvdq=',avgSens_dvdq)
        #print('avgddeldp=',avgSens_ddeldp)
        
    numact=len(act_locs)
    # avgSens_dvdq is list of scalar complex vals
    if parmObj.get_version()==1:  
        q=np.mean(np.absolute(avgSens_dvdq))*numact # take avg across abs value of perf-act pair sensitivities
        p=np.mean(np.absolute(avgSens_ddeldp))*numact
        #print('q=',q) # print when tuning c1 and c2
        #print('(1/q)*c[0]=',(1/q)*c[0])
    else: # in droop case, dont reduce Fq and Fp by number of actuators
        q=np.mean(np.absolute(avgSens_dvdq))*numact # take avg across abs value of perf-act pair sensitivities
        p=np.mean(np.absolute(avgSens_ddeldp))*numact
    try: 
        Fq_ub=(1/q)*c[0]
        Fp_ub=(1/p)*c[1]
    except: 
        print('error: q or p = 0, avgSens are=',avgSens_dvdq)

    return Fq_ub,Fp_ub

              
def detControlMatExistence(parmObj,feeder, act_locs, A, B, indicMat,substation_name,perf_nodes,depths,node_index_map,file_name):
#def detControlMatExistence(feeder, act_locs, perf_nodes,A,B,R,X,indicMat):
    n=int(len(indicMat)/6) # indicMat is 6n x 6n

# Compute good (Fp,Fq) sample space
    R,X=hm.createRXmatrices_3ph(feeder, node_index_map,depths)
    Fq_ub,Fp_ub=computeFParamSpace_v2(parmObj,feeder, act_locs, perf_nodes,R,X,depths,node_index_map,file_name)

    numsamp=15 # temporary, should really do at least 15
    Fq_range=np.linspace(0.0001, Fq_ub, numsamp)
    Fp_range=np.linspace(0.0001, Fp_ub, numsamp)
    print('Fp_range=',Fp_ub,' and Fq_range=',Fq_ub)

# Initialize arrays, will be populated in loops
    feas=False # boolean
    numfeas,myCosts,domeig_mags= (np.empty((0,1)) for i in range(3))
    feasFs,myFbases=(np.empty((0,2)) for i in range(2)) # 0 for dim to concatenate on, 2 for length of vectors being concatenated
    dataFull, data_zero = (np.array([]) for i in range(2)) # for each config
    CLmat=np.empty((6*n,6*n))

    ssError_no_contract,eigs_outside_circle=0,0 # will save which checks on F fail
    
    for Fq in Fq_range:
        for Fp in Fp_range:
            if not(Fq==0 or Fp==0): # skip iteration if either zero
                #print("(Fp,Fq)=",Fp,",",Fq,")")
                if np.isnan(Fp) or np.isnan(Fq):
                    sys.exit('Error: Fq or Fp are NaN')

                F=assignF(parmObj,Fp,Fq,indicMat)
                #print('sizeA=',B.shape)
                CLmat=A-np.dot(B,F) # CLmat=A-BF
                eigs,evecs=LA.eig(CLmat) # closed loop eigenvalues
                eigMags=np.absolute(eigs)

            #     MODIFY: for version 2 (droop), check that evals are within unit circle,
            #     AND vss for vdbc1 and vdbc2 are in 5% range
                if parmObj.get_version()==2: # volt-var and volt-watt control      
                    delV=0.0535 # if dbc causes this change in v, want to ensure that all new steady states are within 0.05
                    dbc_vec1=delV*np.ones((len(CLmat), 1)) 
                    delvss1=np.dot(LA.inv(np.identity(len(CLmat))-CLmat),dbc_vec1) # inv(I-BF)*dbc_vec
    #                 print('vss1=',vss1[:4])
    #                 print('vss1=',vss1[-4:])

                    # identify which nodes have delvss=dbc_vec, and throw those out
                    counter=0
                    delvss1_cleaned=[]
                    for v in delvss1:
                        if delV!=v: # basically only substation node voltages are unaffected by kgains, so we throw those out
                            delvss1_cleaned.append(v)
                        else:
                            counter+=1

                    #print('delvss1 (min,max,median)=(',np.around(np.amin(delvss1_cleaned),5),',',np.around(np.amax(delvss1_cleaned),5),',',np.around(np.median(delvss1_cleaned),5),')')
                    #print('num of immovable 3-ph nodes voltages=',counter/3)
                    bool2=all(item <0.05 for item in delvss1_cleaned) # check that the remaining are all under 0.05
                else:
                    bool2=True # PBC case
                                        
                if all(np.around(eigMags,decimals=6)<=1):
                    # require that all evals=1 have null space full of base
                    # evecs (no generalized evecs)
                    tol=0.0001
                    eval=1
                    num1evals=sum(np.absolute(np.absolute(eigs)-1)<tol) # numel(find(abs(abs(eigs)-1)<tol))   
                    #num1evals=sum(np.absolute(eigs)==1) # numel(find(abs(abs(eigs)-1)<tol))   
                    Y = LA.null_space(CLmat-eval*np.eye(len(CLmat))) #null(CLmat-eval*eye(size(CLmat,1))); % Y is orthonorm basis matrix
                    dimNull=len(Y[0]) # number of cols
                    
                    #print('eigs are in/on unit circle..')
                    #print('num1evals=',num1evals)
                    #print('dimNull=',dimNull)
                    if bool2 and (dimNull==num1evals):                    
                        #print('Found feas F')
                        feasFs=np.append(feasFs,[[Fp, Fq]],axis=0)
                        
                        # find dominant eval
                        mylist = np.absolute(np.absolute(eigs)-1) # find evals not at 1
                        def condition(x): return x > tol # find evals not at 1
                        idx = [i for i, element in enumerate(mylist) if condition(element)] # find evals not at 1
                        if idx: # if nonempty
                            mag_domeig=np.amax(np.absolute(eigs[idx]))
                        else:
                            mag_domeig=eigs[0] # just pick any of them
                        domeig_mags=np.append(domeig_mags,mag_domeig)
                        
                    else: # if F not feas, print in which way it is
                        if not bool2: # save which check failed
                            ssError_no_contract+=1

                        else:
                            eigs_outside_circle+=1

                    
                val=np.sum(eigMags[np.where(eigMags > 1)])
                myCosts=np.append(myCosts,[[val]],axis=0) # temp
                myFbases=np.append(myFbases,[[Fp, Fq]],axis=0) # save all Fs tried
                       
    threshold=1 # your choice, define because if only found 1 feas config too borderline to count as feas
    numfeas=np.append(numfeas,[[len(feasFs)]],axis=0) # number of rows
    #print('feasFs=',np.around(feasFs,3))
    
   # if feas==True:
    if len(feasFs)>=threshold:
        print("Config good!")
        bestF=feasFs[np.argmin(domeig_mags)][:] # the F that results in the most stable dominant eval
        #print("Best F is (Fp Fq)=",bestF) # typically really tiny, not interesting to print
        feas=True
    else:
        feas=False
        print("No F found for config --> bad config")
        print("Unit circle fails=",eigs_outside_circle,', ss_error contraction fails=',ssError_no_contract)

    dataFull=np.concatenate((myFbases,myCosts),axis=1) # [Fp Fq myCost]  
    #print("[Fp,Fq,myCost]=\n",dataFull) # print useful data
    numTried=len(dataFull) # number of rows
    num_act=np.count_nonzero(indicMat)/2
       
    # return feas,feasFs,num_act,numfeas,numTried
    return feas,feasFs,numfeas,numTried,num_act