function H = heavisideby(x, k)
%HEAVISIDEBY Summary of this function goes here
%   Detailed explanation goes here

H = 1/(1 + exp(-2*k*x));

end

