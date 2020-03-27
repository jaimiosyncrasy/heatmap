function [P12,solns] = makePVcurve(Sweep_lb,Sweep_ub,Sbase,Vbase,R12,X12,V1)
% varying P, so plot PV and Pdel curve

numPts=20; % number of points to comprise curves
P12=linspace(Sweep_lb,Sweep_ub,numPts); % not pu
Q12pu=tan(acos(0.9)) % center Q so that when P12=Sbase pow factor=0.9
Q12=Q12pu*Sbase; % set const, not pu

    trueV2=[]; trueDel2=[];
    for i=1:length(P12)
        % Solve nonlin eqns with fwd bwd sweep
        [a,b] = solveFwdBwdSweep_2bus(R12,X12,V1,P12(i),Q12);
        trueV2=[trueV2 a];
        trueDel2=[trueDel2 b]; % degrees
        % Solve lzn PF eqns
        syms V2 delta2
        eqn1= (abs(V1))^2-V2^2==2*R12*P12(i)+2*X12*Q12;
        eqn2= angle(V1)-delta2==(X12*P12(i)-R12*Q12)/(abs(V1)*V2);
        sol=solve([eqn1, eqn2],[V2 delta2])
        lznV2(i)=eval(sol.V2(1)); % degrees
        lznDel2(i)=(180/pi)*eval(sol.delta2(1));
    end

    %% plot lzn curve with true curve
    figure;  hold on;
    plot(P12/Sbase,lznV2/Vbase,'r-',P12/Sbase,trueV2/Vbase,'b-','LineWidth',2); xlabel('P12, kW'); ylabel('V2, pu'); legend('linearization','true');
    title('True P-V Curve and Linearization Curve'); ylabel('Voltage, pu');
    
    figure; hold on;
    plot(P12/Sbase,lznDel2,'r-',P12/Sbase,trueDel2,'b-','LineWidth',2); xlabel('P12, kW');  ylabel('Delta2, degrees'); legend('linearization','true');
    title('True P-Del Curve and Linearization Curve'); ylabel('Angle, degrees');

    % make structure so that solns holds all 4 vectors
    solns.trueV2=trueV2; solns.trueDel2=trueDel2; solns.lznV2=lznV2; solns.lznDel2=lznDel2;
end

