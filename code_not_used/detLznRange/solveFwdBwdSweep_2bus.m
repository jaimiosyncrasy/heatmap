function [V2,del2] = solveFwdBwdSweep_2bus(R12,X12,V1,P2,Q2)
% Solves power flow for 2-bus grid using fwd backward sweep iteration
% method
% Based on code from text pg 346, Ex 10.1, 3 bus (one slack, two PQ load)

%% Initialization
disp('~~~~~~~ Starting FBS Method for Solving PF');

% Givens: z12, V1, P2, Q2
S2=P2+j*Q2; % per unit
z12=R12+j*X12; % per unit
Vs =V1; % per unit

% Init Cond
Vnom = Vs; % to check convergence
tol = 0.0001;
k = 1; % rmb matlab indices start at 1
V1(k) = 0;
V2(k) = 0; 
Vconv{k}=[0 0];

%% First Iteration
k = k+1; % new iteration
% Fwd Sweep
V1(k) = Vs;
V2(k) = Vs;

% check convergence:
Vconv{k} = [abs((abs(V1(k))-abs(V1(k-1))))/Vnom ...
abs((abs(V2(k))-abs(V2(k-1))))/Vnom];

% Backward sweep
I12 = (S2/V2(k))';

%% Iterative Part

while(any(Vconv{k}>=tol)) % break when all nodes less than tol, means continue while any of the nodes are greater 
    k = k+1;  % new iteration
    % Fwd sweep
    V1(k)=V1(k-1); % same as prev iter
    V2(k)=Vs-z12*I12;

    % check convergence:
    Vconv{k} = [abs((abs(V1(k))-abs(V1(k-1))))/Vnom ...
    abs((abs(V2(k))-abs(V2(k-1))))/Vnom];
    % Vconv{:} uncomment when debugging
    
    % Backward sweep
    I12 = (S2/V2(k))';
    
    if length(Vconv)>30
        disp('Didnt converge');
        break; % break out of loop if iter too many times
    end

end

%% Output Results
disp('~~~~~~~ PF Results: ');
Vsoln = [V1(end) V2(end)]  % /Vs to put into pu
convergedIfZero=Vconv{end}
numIter = length(Vconv)-1 % -1 because Vconv(1) initialized at zero
disp('~~~~~~~ Finished FBS Method for Solving PF');

%% Polar to rect conversion for testing/probing
foo = Vsoln;
Vsoln_polarCoor = [[abs(foo)].' [angle(foo)*180/pi].'] % Vpu, deg
V2=abs(Vsoln(2)); % 
del2=angle(Vsoln(2))*180/pi; % degrees
end

