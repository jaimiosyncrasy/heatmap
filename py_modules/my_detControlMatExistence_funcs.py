import importlib
import setup_nx # your own module, setup.nx.py
import numpy as np
import scipy.linalg as LA
import math
import statistics as st
import cmath
import matplotlib.pyplot as plt 
import itertools
from operator import add
importlib.reload(setup_nx)
from setup_nx import *
from graphviz import Source, render
from sympy import * # includes Matrix object and solve function

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

def assignF_v2(f_ub_table,n): # algo similar to updateStateSpace
    # create one F matrix from parm space
    # indicMat has diagonal 3x3 blocks where an APNP is
    F=np.zeros((6*n,6*n)) # make numpy matrix that will hold values not symbolic var
    candFset,candFset_percentUB=np.zeros((1,len(f_ub_table))),np.zeros((1,len(f_ub_table)))
        
    for k in range(0,len(f_ub_table)): # for each actuator
        # fubtable has format [APNP_idx indicMat_row indicMat_col flb fub]
        fub=f_ub_table[k][4] # 3 for 4th col, which has upper bound
        #print('fub=',fub)
        i=int(f_ub_table[k][1]) # i from nonzero Fij
        j=int(f_ub_table[k][2]) # j from nonzero Fij
        sigma=1 # arbitrary variance, may need to adjust
        #probs=[0.05, 0.1, 0.1, 0.25, 0.45, 0.05, 0.0, 0.0, 0.0, 0.0] # probability distr, elements must sum to 1
        probs=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1] # probability distr, elements must sum to 1
        #probs=[0, 0, 0, 0, 0.05, 0.1, 0.1, 0.25, 0.45, 0.05] # probability distr, elements must sum to 1
        fval=np.random.choice(np.linspace(0,fub,len(probs)), p=probs) # take sample from backwards weibul, 45% chance it's 90% of fub
        #print('fval=',fval)
        candFset[0][k]=fval # save fvals into vec
        candFset_percentUB[0][k]=fval/fub 
        # assume all blocks are 3ph for now, will correct later in this funcc
        #print('(i,j)=',i,' and ',j)
        F[i,j]=fval

    return F,candFset,candFset_percentUB
    
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
    if '13NFbalanced' in file_name:
        c=np.array([0.412,0.857]) # (q,p) tuned for 13NF based on data from feas configs
    elif '123NF' in file_name:
        #c=np.array([0.3,0.45])  # setting on 9/2/20
        #c=np.array([0.5,0.7]) # place_max_coloc getting not-great results from this
        if parmObj.get_version()==1:  
            c=np.array([0.3,0.45])  # setting back to this on 2/28/21
        else: # droop case
            c=np.array([0.7,1.1])  # setting this on 3/7/21
    elif 'PL0001' in file_name: 
        c=np.array([0.32,0.65])
    elif 'oaklandJ' in file_name:
        c=np.array([0.3,0.45])
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
    zpath_var1=np.var(np.absolute(avgSens_dvdq)) # if z to sub varies a lot per act, more stable so dec p and q
    zpath_var2=np.var(np.absolute(avgSens_ddeldp))
    print('zpath_multiplier1=',np.around(zpath_var1,3))
        
    m=lambda var,nact: 1+(nact-1)*np.exp(-120*var) # lambda func, for high variance m=1, for small variance m=numact
    q=np.mean(np.absolute(avgSens_dvdq))*m(numact,zpath_var1) # take avg across abs value of perf-act pair sensitivities
    p=np.mean(np.absolute(avgSens_ddeldp))*m(numact,zpath_var1)
    #print('q=',q) # print when tuning c1 and c2
    #print('(1/q)*c[0]=',(1/q)*c[0])

    try: 
        # added a floor kgain of (0.05,0.1) due to realism of kgain settings. Result is all controllers will have a limit to # DERs they can place
        Fq_ub=max((1/q)*c[0],0.05) # lower bound on kgains to sample
        Fp_ub=max((1/p)*c[1],0.1)
    except: 
        print('error: q or p = 0, avgSens are=',avgSens_dvdq)

    return Fq_ub,Fp_ub

def computeFParamSpace_v3(B_nonmat,indicMat,indicMat_table):
    B = Matrix(B_nonmat)# convert to matrix opject so can do B.multiply(Fsymb) later
    n=int(len(indicMat)/6) # indicMat always has 6n rows
    row,col=np.nonzero(indicMat) # find indices of nonzero f ele
    #print('row=',row,' and col=',col)
    # if apply gdisc cond to every nonzero ele in F
    f_loc=indicMat_table[:,1:3] # get last 2 cols
    #np.concatenate((np.transpose([row]),np.transpose([col])),axis=1)
    # if only apply gdic cond to one f ele per 3ph block: 
    #f_loc=np.concatenate((np.transpose([row[::3]]),np.transpose([col[::3]])),axis=1) # only get index of first 1 in diagonal 3x3 block     
   
    numact=len(np.unique(f_loc[:,0])) # if one 3ph actuator, numact=3
    numperf=len(np.unique(f_loc[:,1]))      
    num_fele=len(f_loc) 
    print('num_fele=',num_fele)

    # initialize c,r, and bigBox_range
    c=np.zeros((6*n,1))
    c=Matrix(c) # convert to matrix obj so can have symbolic ele
    r=np.zeros((6*n,1))
    r=Matrix(r) 

    # main function call:
    bigBox_range=np.around(det_Gdisc_Franges(indicMat_table, n, B, c, r),2) # round to two decimal places

    return bigBox_range
    

# helper func in computeFParamSpace_v3
def multBF(B,i,j,n):
    # returns B*F faster than if did B.multiply(Fsymb)
    # uses how the only nonzero ele of F is at Fij, so only one col of BF is nonzero
    Bcol=B[:,i]
    var('f') # declare symbolic
    BF=zeros(6*n) # faster than creating array and converting to matrix object
    BF[:,j]=Bcol*f
    return BF

# helper func in computeFParamSpace_v3, calls det_Gdisc and solve_Gdisc_conds
def det_Gdisc_Franges(indicMat_table, n, B, c, r):
# for each row of matrix BF, determine Gdisc conditions in terms of radius and center
# reduces B and F matrix before applying gdisc conditions
# solve gdisc conditions to compute range on F gains, store in bigBox
    f_loc=indicMat_table[:,1:3] # get last 2 cols
    num_fele=len(f_loc) 
    var('f') # declare symbolic
 
    bigBox_range= np.concatenate((indicMat_table, np.zeros((len(indicMat_table),2))), axis=1) # add 2 columns to end
    
    for k in range(num_fele): # apply gdisc conditions for each nonzero f ele, each iteration updates row of bigBox_range
        print('working with ',k,'_th nonzero F ele')
        Fsymb=zeros(6*n) # faster than creating array and converting to matrix object
        #num_col_indicMat=len(indicMat[0]) # number of cols in indicMat
        injNode=f_loc[k][0]
        sensNode=f_loc[k][1]
        if injNode<3*n and sensNode>3*n: 
            raise Exception("cant choose f_loc in upper right of F, f_loc=["+str(injNode)+" "+str(sensNode)+"]")
        else:
            Fsymb[injNode,sensNode]=f
            #c[k]=myBF[sensNode,sensNode] # f_ij has gdisc center at BF_{jj} because BF's jth col is nonzero
        myBF=multBF(B,injNode,sensNode,n)
        nnz_idx=0 # initialize index of nonzero center, should only have one because look at one f ele at a time

        c,r,nnz_idx=det_Gdisc(c,r,myBF,nnz_idx,n) 
        print('solving Gdisc conditions...')
        bigBox_range=solve_Gdisc_conds(c,r,nnz_idx,k,bigBox_range) # updates bigBox range with a new range
    
    # print('bigBox_range=',bigBox_range)
    return bigBox_range
    
# helper func in computeFParamSpace_v3
def det_Gdisc(c,r,myBF,nnz_idx,n):
# for each row of matrix BF, determine radius and center of gdisk
# BF should be full size, not reduced; if reduced, rdc_BF[i,i] has no meaning
    for i in range(6*n): # one gdisc for each row of BF (pxp)
        c[i]=myBF[i,i] # index a Matrix object with [row,col] whereas np array with [row][col]
        if c[i]!=0 and nnz_idx==0: 
            nnz_idx=i # should only ever be one because there is one f ele
        sum=0
        
        for j in range(6*n): # for each col of F, F is axp
            #sum=sum+myBF_noD[i,j]
            if i!=j:
                sum=sum+np.absolute(myBF[i,j]) # row sum
            #print(sum)
        r[i]=sum
    return c,r,nnz_idx

# helper func in computeFParamSpace_v3
def solve_Gdisc_conds(c,r,nnz_idx,k,bigBox_range):
    # must pass in bixBox_range so can update it across calls of this func
    # modify indicMat_table by adding 2 columns to end, for [flb fub] computed from gdisc cond
                 
    #print('r=',r)
    #print('c=',c)
    cond1=c[nnz_idx]+r[nnz_idx]-2 # =0, c+r<2
    print('cond1=',cond1)
    cond2=c[nnz_idx]-r[nnz_idx] # =0, c-r>0
    #print('cond2=',cond2)
    cond3=c[nnz_idx] # =0, c>0
    #print('cond3=',cond3)

#     sol1=solve(cond1,dict=True) # eqn is set to zero and solved
#     sol2=solve(cond2,dict=True) 
#     if sol1 and sol2: # if both give solns
#         if sol1[0][f]<sol2[0][f]:
#             bigBox_range[k][3:5]=[sol1[0][f], sol2[0][f]] # row of bigBox_range
#         else:
#             bigBox_range[k][3:5]=[sol2[0][f], sol1[0][f]]
#         #print('range=',bigBox_range[k][:])
#         bigBox_range[k][0:3]=indicMat_table[k]
#     else:
#         print('No solution to conditions')
     
# only solve first condition because second condition is always cf>0 --> f>0
    sol1=solve(cond1,dict=True) # eqn is set to zero and solved
    if sol1: # if both give solns
        bigBox_range[k,3:5]=[0, sol1[0][f]] # last 2 cols
    else:
        print('No solution to conditions') 
        
    return bigBox_range

#-------------------------------------------------------------------------------------------------------------
    
def eval_Fmat(parmObj,CLmat,domeig_mags): # items we populate for each F tried
    ssError_no_contract,eigs_outside_circle,Fstable=0,0,0 # saves which checks on F fail

    eigs,evecs=LA.eig(CLmat) # closed loop eigenvalues
    eigMags=np.absolute(eigs)

#     For version 2 (droop), check that evals are within unit circle,
#     AND vss for vdbc1 and vdbc2 are in 5% range
    if parmObj.get_version()==2: # volt-var and volt-watt control      
      #  bool2=True # temporarily removing from droop
        
        delV=0.0535 # if dbc causes this change in v, want to ensure that all new steady states are within 0.05
        dbc_vec1=delV*np.ones((len(CLmat), 1)) 
        delvss1=np.dot(LA.inv(np.identity(len(CLmat))-CLmat),dbc_vec1) # inv(I-BF)*dbc_vec

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
            Fstable=1 # means F is stabilizing 

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
                ssError_no_contract=1

            else: # nullity condition doesnt hold
                eigs_outside_circle=1

    else: # outside unit circle
        eigs_outside_circle=1

    val=np.sum(eigMags[np.where(eigMags > 1)])
    return val,ssError_no_contract,eigs_outside_circle,domeig_mags,Fstable 
    
    
def detControlMatExistence(parmObj, feeder, A, B, indicMat,indicMat_table,act_locs,perf_nodes,node_index_map,depths,file_name):
#def detControlMatExistence(feeder, act_locs, perf_nodes,A,B,R,X,indicMat):
    n=int(len(indicMat)/6) # indicMat is 6n x 6n
    print('indicMat size is ',n)

    #sample_way='heuristic' # you choose
    sample_way='non-heuristic' # you choose
    
    # Initialize arrays, will be populated in loops
    feas=False # boolean
    numfeas,myCosts,domeig_mags= (np.empty((0,1)) for i in range(3))
    dataFull, data_zero = (np.array([]) for i in range(2)) # # [Fp Fq myCost] for each config
    CLmat=np.empty((6*n,6*n))
    eigs_outside_circle,ssError_no_contract=0,0
    print('evaluating kgains sampled from parm space...')
    numsamp=225 # temporary, should really do at least 15

    #-------------------------------------------------
    if sample_way=='heuristic':
        # Heuristically compute sample space
        R,X=hm.createRXmatrices_3ph(feeder, node_index_map,depths,file_name) # need for computing parm space ranges
        Fq_ub,Fp_ub=computeFParamSpace_v2(parmObj,feeder, act_locs, perf_nodes,R,X,depths,node_index_map,file_name)
        print('Fp_range=',Fp_ub,' and Fq_range=',Fq_ub)
        assert(math.sqrt(numsamp) ** 2 == numsamp) # checks that numsamp is a perfect square
        Fq_range=np.linspace(0.0001, Fq_ub, round(sqrt(numsamp)))
        Fp_range=np.linspace(0.0001, Fp_ub, round(sqrt(numsamp)))  
        
        feasFs=np.array([], dtype=np.int64).reshape(0,2)
        myFbases=np.array([], dtype=np.int64).reshape(0,2)
        
        # heuristic iter
        for Fq in Fq_range:
            for Fp in Fp_range:
                if not(Fq==0 or Fp==0): # skip iteration if either zero
                    #print("(Fp,Fq)=",Fp,",",Fq,")")
                    if np.isnan(Fp) or np.isnan(Fq):
                        sys.exit('Error: Fq or Fp are NaN')

                    F=assignF(parmObj,Fp,Fq,indicMat)
                    #print('sizeA=',B.shape)
                    CLmat=A-np.dot(B,F) # CLmat=A-BF
                    candFset=np.zeros((1,2))
                    candFset[0][0],candFset[0][1]=Fq,Fp # save fvals into vec
                    val,ssErr_bool,ecirc_bool,domeig_mags,Fstable=eval_Fmat(parmObj,CLmat,domeig_mags)

                    if Fstable:
                        print('found stabilizing F!')
                        feasFs=np.append(feasFs,candFset,axis=0) # if Fstable, add Fset to new row of feasFs
                    eigs_outside_circle+=ecirc_bool
                    ssError_no_contract+=ssErr_bool
                    myCosts=np.append(myCosts,[[val]],axis=0) # temp
                    myFbases=np.append(myFbases,candFset,axis=0) # save all Fs tried, not just feas ones
   
    #-------------------------------------------------
    elif sample_way=='non-heuristic':

        # Non-heuristic way to compute the sample space
        f_ub_table=computeFParamSpace_v3(B, indicMat,indicMat_table) # format is [Frow Fcol f_lb f_ub]
        str_arr=np.array([], dtype=np.int64).reshape(0,2)
        print('ex=',[act_locs[int(0)], perf_nodes[int(0)]])
        for idx in f_ub_table[:,0]: # pull APNP idx from first col of f_ub_table
            str_arr=np.append(str_arr,[np.array([act_locs[int(idx)], perf_nodes[int(idx)]])],axis=0) # adds 2x1 col for each APNP pair
        print('parm space table=\n',np.append(str_arr,f_ub_table[:,1:],axis=1),' << formated as [act perf row col lb ub]')
        feasFs=np.array([], dtype=np.int64).reshape(0,len(f_ub_table)) # number of cols will be = number of nonzero F ele
        myFbases=np.array([], dtype=np.int64).reshape(0,len(f_ub_table))

            
        # non-heuristic iter
        for k in range(numsamp):
            F,candFset,candFset_percentUB=assignF_v2(f_ub_table,n) # design F matrix
            CLmat=A-np.dot(B,F) # CLmat=A-BF
            val,ssErr_bool,ecirc_bool,domeig_mags,Fstable=eval_Fmat(parmObj,CLmat,domeig_mags)
            print('cand F set=',candFset)
            if Fstable:
                print('found stabilizing F!')
                feasFs=np.append(feasFs,candFset,axis=0) # if Fstable, add Fset to new row of feasFs
            eigs_outside_circle+=ecirc_bool
            ssError_no_contract+=ssErr_bool
            myCosts=np.append(myCosts,[[val]],axis=0) # temp
            myFbases=np.append(myFbases,candFset,axis=0) # save all Fs tried, not just feas ones

        #-------------------------------------------------
    else:
        raise Exception("unrecognized string value for sample_way")
        
# [for both heuristic and non-heuristic] gather result
    threshold=(0.1)*(numsamp) # your choice, defines when node color is yellow vs. blue
    numfeas=np.append(numfeas,[[len(feasFs)]],axis=0) # number of rows
    #print('feasFs=',np.around(feasFs,3))
    
   # if feas==True:
    if len(feasFs)>=threshold:
        print("Config good!")
        print('feasF shape=',feasFs.shape)
        bestF=feasFs[np.argmin(domeig_mags)][:] # the (Fp,Fq) that results in the most stable dominant eval
        #print("Best F is (Fp Fq)=",bestF) # typically really tiny, not interesting to print
        feas=True
    else:
        bestF=float("NaN")
        feas=False
        print("No F found for config --> bad config")
        print("Unit circle fails=",eigs_outside_circle,', ss_error contraction fails=',ssError_no_contract)

    dataFull=np.concatenate((myFbases,myCosts),axis=1) # [Fp Fq myCost]  for all tried Fs
    #print("[Fp,Fq,myCost]=\n",dataFull) # print useful data
    numTried=len(myCosts) # number of rows
    num_act=np.count_nonzero(indicMat)/2
       
    # return feas,feasFs,num_act,numfeas,numTried
    return feas,feasFs,numfeas,numTried,num_act,bestF,indicMat

def eval_config_F(parmObj,bestFparm,indicMat,feeder,depths,node_index_map):
                    # NEED FIX
    # takes in particular F mat +indicMat and evaluates whether (A-BF) is stable; if so feas=True
    Fp=bestFparm[0]
    Fq=bestFparm[1]
    F=assignF_v2(parmObj,Fp,Fq,indicMat)
    
    n=int(len(F)/6) # indicMat and F is 6n x 6n
    A, B = hm.setupStateSpace(parmObj,feeder,n,node_index_map,depths)

    # simpler version of detControlMatExistence:
    # Initialize arrays, will be populated by eval_Fmat function
    domeig_mags= (np.empty((0,1)) for i in range(1))
    feasFs=np.array([], dtype=np.int64).reshape(0,2)
    ssError_no_contract,eigs_outside_circle=0,0 # will save which checks on F fail

    CLmat=np.empty((6*n,6*n))
    CLmat=A-np.dot(B,F) # CLmat=A-BF
    val,ssError_no_contract,eigs_outside_circle,domeig_mags,feasFs=eval_Fmat(parmObj,CLmat,candFset,
                                                                             ssError_no_contract,eigs_outside_circle,domeig_mags,feasFs)
    # feasFs=0 if F unstable, feasFs=1 if F is stable
    print('len feasFs=',len(feasFs))
    if len(feasFs)>=1:
        feas=True
    else:
        feas=False

    print('Actuator configuration is feasible') if feas else print('Actuator configuration is not feasible')
    return feas