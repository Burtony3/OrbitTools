function [H, stats] = M2H(M, e, units)

% INPUT HANDLING
if nargin < 3; units = 'degrees'; end
units = char(units);
if strcmpi(units, 'd')
    units = 180/pi;
else
    units = 1;
end
M = M / units;

% ERROR
error = 1.e-8;
 
% INITIALIZING:
H = sign(M)*min(abs(M), 10);
 
% LOOPING
ratio = Inf;
iters = 0;
HList = H;
while abs(ratio) > error
    iters = iters+1;
    ratio = (e*sinh(H) - H - M)/(e*cosh(H) - 1); % Newton Step
    H = H - ratio;
    HList(end+1) = H*units; %#ok<AGROW>
end

% OUTPUTTING
stats.Iterations = iters;
stats.H = HList;
H = H*units;

end