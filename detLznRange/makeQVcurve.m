function [Q12,solns] = makeQVcurve(Sweep_lb,Sweep_ub,Sbase,Vbase,R12,X12,V1)
% varying Q, so plot QV and Qdel curve

numPts=20; % number of points to comprise curves
Q12=linspace(Sweep_lb,Sweep_ub,numPts); % not pu
P12pu=tan(acos(0.9)) % center Q so that when P12=Sbase pow factor=0.9
P12=P12pu*Sbase; % set const, not pu

    trueV2=[]; trueDel2=[];
    for i=1:length(Q12)
        % Solve nonlin eqns with fwd bwd sweep
        [a,b] = solveFwdBwdSweep_2bus(R12,X12,V1,P12,Q12(i));
        trueV2=[trueV2 a];
        trueDel2=[trueDel2 b]; % degrees
        % Solve lzn PF eqns
        syms V2 delta2
        eqn1= (abs(V1))^2-V2^2==2*R12*P12+2*X12*Q12(i);
        eqn2= angle(V1)-delta2==(X12*P12-R12*Q12(i))/(abs(V1)*V2);
        sol=solve([eqn1, eqn2],[V2 delta2])
       lznV2(i)=eval(sol.V2(1)); % degrees
        lznDel2(i)=(180/pi)*eval(sol.delta2(1));
    end

    %% plot lzn curve with true curve
    figure; hold on;
    plot(Q12/Sbase,lznV2/Vbase,'r-',Q12/Sbase,trueV2/Vbase,'b-','LineWidth',2); xlabel('Q12, kW'); ylabel('V2, pu'); legend('linearization','true');
    title('True Q-V Curve and Linearization Curve'); ylabel('Voltage, pu');
    
    figure; hold on;
    plot(Q12/Sbase,lznDel2,'r-',Q12/Sbase,trueDel2,'b-','LineWidth',2); xlabel('Q12, kW'); ylabel('Delta2, degrees'); legend('linearization','true');
    title('True Q-Del Curve and Linearization Curve');

    % make structure so that solns holds all 4 vectors
    solns.trueV2=trueV2; solns.trueDel2=trueDel2; solns.lznV2=lznV2; solns.lznDel2=lznDel2;
end

