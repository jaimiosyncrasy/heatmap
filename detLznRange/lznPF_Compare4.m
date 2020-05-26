clc; close all; clear all;
% All values are in pu
% From 13NF PF results, substation is delivering this much to the feeder: 1250 kW, 680kVAR
% on 13NF, typical loads are 100*10^3 W, 50*10^3 VAR
% a typical line is 0.07+0.21j ohms

% Inductive line: with pos+j*pos impedance, P+jQ is pos+pos*j too, though del2=neg

%%
%1) V1=1, del1=0, slack bus
Vbase=4160/sqrt(3); % medium voltage, 13NF primary side voltage, not in per unit
V1=1*Vbase; % slack, not pu (4160 volt feeder)

%2) z12 is the sum of impedanced between the substation and the performance node (not actuator node)
z12=(0.15+0.5*j); % inductive, ohms, subst to  bottom part of 13NF, no pu
R12=real(z12); X12=imag(z12);

%3) V2 and del2 will computed from solving power flow on the 2-bus equivalent system

%4) P12 and Q12 is varied in some range. For each value in this range we solve power flow.
% So we need to define ranged for these "sweeps". Then we should define an Sbase value, arbitrarily taken as the mean of the sweep bounds.
% See overleaf for details

%% Define power curve sweep range
% Without the full distribution grid network, in this code we use some hacky way to come up with sweep values
% You can start with this implementation, then go back and use the feeder object to come up with Sbase according
% to procedure in overleaf StateSpace pdf, section 7

V2=(0.98*cos(-5*pi/180)+j*0.98*sin(-5*pi/180))*Vbase
Iest=(V1-V2)/z12
Ibase=round(mean([abs(real(Iest)),abs(imag(Iest))]),-1)
Zbase=Vbase/Ibase
Sbase=Vbase^2/Zbase % Watts, use to define sweep of P & Q

Sweep_lb=0.5*Sbase;
Sweep_ub=1.5*Sbase;

%% Make PV curves 
% vary P, plot P-V and P-del curve
[pvals,solns1] = makePVcurve(Sweep_lb,Sweep_ub,Sbase,Vbase,R12,X12,V1);

% vary Q, plot Q-V and Q-del curve
[qvals,solns2] = makeQVcurve(Sweep_lb,Sweep_ub,Sbase,Vbase,R12,X12,V1);

%% Run opt problem that solves for lzn itvl
[errVmax1,slope_vp]=computeLznItvl2(pvals/Sbase,solns1.lznV2/Vbase,solns1.trueV2/Vbase,1); % (x,fx_lzn,fx_true)
[errDelmax1,slope_delp]=computeLznItvl2(pvals/Sbase,solns1.lznDel2,solns1.trueDel2,2); % (x,fx_lzn,fx_true)

[errVmax2,slope_vq]=computeLznItvl2(qvals/Sbase,solns2.lznV2/Vbase,solns2.trueV2/Vbase,3); % (x,fx_lzn,fx_true)
[errDelmax2,slope_delq]=computeLznItvl2(qvals/Sbase,solns2.lznDel2,solns2.trueDel2,4); % (x,fx_lzn,fx_true)

%% Outputs: for each curve (4 types), 

% find max true-lzn error, 
[errVmax1, errVmax2, errDelmax1, errDelmax2] % 1 is when vary P, 2 is when vary Q
% errVmax is in Vpu, errDelmax is in deg

% find max/min/mean slope across P/Qsweep for true curve
[slope_vp, slope_delp, slope_vq, slope_delq]
% slopes are in Vpu/Spu or Degrees/Spu

%%
% Full Computation of linearization procedure:
% 1) Input: network (impedances, loads), actuator location, perf node
% location
% 2) Call method computeZeff
% 3) Set v1=1, delta1=0
% 5) Determine P & Q sweep ranges, and Sbase
% ------- For PV,Pdel,QV,Qdel...
% 6) solve power flow multiple times for 2-bus equiv network
% 7) Draw lzn and true curves
% 8) Call method computeItvl
% -----------------------------
% 9) Set u capacity and determine dbc bounds
% Output: lzn intervals and dbc bounds wrt given inputs