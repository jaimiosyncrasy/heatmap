%% 6.3.20 Compute Line Losses for first 8-node network
% inputs: set of P at TOD, set of Q at TOD, feeder network object, act locs (ax1 vector),
    % lb_mode (boolean), actuator capacity (ax1 vector per unit) 
    % take nominal actuator capacity as 600kVAR so get pu by dividing that by
    % Sbase
    % note: dont need performance node, but do need actuator nodes (what I said
    % previously was wrong)
% outputs: Peq,Qeq
% where fit in with detLznRange? called in computePQsweep

clc; clear all; close all;
%--------------------------------------------------
% Setup variables ---------------------------------
n=7;
P=0.5*ones(1,n); % real power at all nodes except slack bus
Q=0.2*ones(1,n); % reactive power
S=P+j*Q;
z=(0.02+0.05j)*ones(n,n); % impedances, get from feeder network object
lb_mode=true;

% guess Vest as a value proportional to distance from substaion.
zmax=10; % longest impedance to substation, in reality you'd compute it using impedance path function
zpath4=7; zpath6=10; zpath7=9; % impedance path to each edge node of nwk
Vest(4)=1-0.2*(zpath4/zmax); % per unit
Vest(6)=1-0.2*(zpath6/zmax); 
Vest(7)=1-0.2*(zpath7/zmax); 

% When looking at computing Plb or Qlb, actuator worst case scenario is
% netative the capacity value. WHen looking at computing
computing_lb=true; % pass computing_lb argument in, so that you can say soemthing like 

if lb_mode==true 
    Sact(1)=-1; % per unit, actuator capacity value, at node 4
    Sact(2)=-1; % capacity value, at node 5
else
    Sact(1)=1;
    Sact(2)=1;
end
% For example comput impedance between node and subst, and let zmax be
% longest path to substation
V=1-z/zmax
I=zeros(n,n); Sloss=zeros(n,n); % will populate 

% -----------------------------------------
% Compute line losses below
% -----------------------------------------
% Line loss of line 34, edge line
I(3,4)=conj((S(4)+Sact(1))/Vest(4)); % take complex conjugate, Sact1 is at node 4
Sloss(3,4)=I(3,4)^2*z(3,4);

% Line loss of line23, non-edge line
V(3)=V(4)+I(3,4)*z(3,4);
I(2,3)=conj(S(3)/V(3))+I(3,4); % I30+I(3,4)
Sloss(2,3)=I(2,3)^2*z(2,3);

% Notice node 2 is a fork, so we've finished finding losses in branch 2-4

% Line loss of line 56, edge line
I(5,6)=conj(S(6)/Vest(6)); % take complex conjugate
Sloss(5,6)=I(5,6)^2*z(5,6);

% Notice node 5 is a fork, so we've finished finding losses in branch 5-6

% Line loss of line 57, edge line
I(5,7)=conj(S(7)/Vest(7)); % take complex conjugate
Sloss(5,7)=I(5,7)^2*z(5,7);

% Notice node 5 is a fork, so we've finished finding losses in branch 5-7

% Notice among the forks, fork 5 is further downstream than fork 2 so
% resolve fork at node 5 first

% Line loss of line25, line connected to fork at node 5
Vfork1=V(6)+I(5,6)*z(5,6);
Vfork2=V(7)+I(5,7)*z(5,7);
checkVclose=(Vfork1-Vfork2)/Vfork1<0.1; % check not more than 10% different
if checkVclose==false
    error('choose different edge V'); % dont need to implement this error call in python, just manually check that they are close
end
V(5)=mean([Vfork1 Vfork2]) % average of all Vfork estimates, one for each branch of the fork
I(2,5)=I(5,6)+I(5,7)+conj((S(5)+Sact(2))/V(5)); % Sact2 is at node 5
Sloss(2,3)=I(2,5)^2*z(2,5);

% Line loss of line12, line connected to fork at node 2
Vfork1=V(5)+I(2,5)*z(2,5);
Vfork2=V(3)+I(2,3)*z(2,3);
checkVclose=(Vfork1-Vfork2)/Vfork1<0.1; % check not more than 10% different
if checkVclose==false
    error('choose different edge V');
end
V(2)=mean([Vfork1 Vfork2]) % average of all Vfork estimates, one for each branch of the fork
I(1,2)=I(2,3)+I(2,5)+conj(S(2)/V(2));
Sloss(1,2)=I(1,2)^2*z(1,2);

% Line loss of line01, non-edge line
V(1)=V(2)+I(1,2)*z(1,2);
I(2,3)=conj(S(3)/V(3))+I(3,4); % I30+I(3,4)
Sloss(2,3)=I(2,3)^2*z(2,3);

% Now sum everything to get Seq
idxDown=[3 4 5 6 7]; % indices of nodes downstream from perf_node
Seq=sum(sum(Sloss))+sum(S)+sum(Sact);
Peq=real(Seq)
Qeq=imag(Seq)