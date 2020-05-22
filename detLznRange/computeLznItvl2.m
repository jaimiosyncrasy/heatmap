function [err_max,slopeInfo]=computeLznItvl(x,fx_lzn,fx_true,figNum)

% Compute lzn itvl
r=abs(fx_lzn-fx_true);
err_max=max(r);

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
