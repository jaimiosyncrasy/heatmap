function [xlb,xub,fx_max,slopeInfo]=computeLznItvl(x,fx_lzn,fx_true,figNum)

% choose so that curves align, somewhat arbitrary but directly affects opt soln
e=0.15*(max(fx_true)-min(fx_true)) % we say lzn holds if true and lzn are within 15%

% Compute lzn itvl
r=abs(fx_lzn-fx_true);
temp=find(r<e); % indices
lb=min(temp); ub=max(temp);
mean_idx=round(length(x)/2)

% Plot lzn itvls
figure(figNum)
plot([x(lb) x(lb)],[fx_lzn(lb) fx_true(lb)],'k-o',[x(ub) x(ub)],[fx_lzn(ub) fx_true(ub)],'k-o',[x(mean_idx) x(mean_idx)],[fx_lzn(mean_idx) fx_true(mean_idx)],'m-'); 

xlb=x(lb); xub=x(ub);
fx_max=max(fx_true);

%Compute slopes everywhere on true curve
% if were looping, want (fx_true(i-1)-fx_true(i+1))/(x(i-1)-x(i+1))
r=([fx_true 0 0]-[0 0 fx_true])./([x 0 0]-[0 0 x]); % all slopes except edges
s=([fx_true 0]-[0 fx_true])./([x 0]-[0 x]); % edge slopes
slopes_true=[s(2) r(3:end-2) s(end-1)];
slopeInfo=[min(slopes_true) max(slopes_true) mean(slopes_true)];
end

% Test code for slope computation:
    % clc;
    % fx_true=[1 1 2 2 3 3 5 5 6 6 3 3 1 1];
    % x=1:length(fx_true);
    % figure; plot(x,fx_true,'r-o'); hold on;
    % 
    % r=([fx_true 0 0]-[0 0 fx_true])./([x 0 0]-[0 0 x]); % all slopes except edges
    % s=([fx_true 0]-[0 fx_true])./([x 0]-[0 x]); % edge slopes
    % slopes_true=[s(2) r(3:end-2) s(end-1)]
    % axis equal
    % plot(x,slopes_true,'bo');
