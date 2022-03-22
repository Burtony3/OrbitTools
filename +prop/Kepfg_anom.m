function [rf, vf] = Kepfg_anom(r0, v0, dnu, mu, units)
if nargin < 5; units = 'degrees'; end

% ALIGNING UNITS
units = char(units);
if strcmpi(units(1), 'd')
    units = 180/pi;
else
    units = 1;
end
dnu = dnu/units;

% PREALLOCATING
[rf, vf] = deal(zeros(length(dnu), 3));


% PULLING CONSTANTS
R0 = norm(r0);
% V0 = norm(v0);
h = norm(cross(r0, v0));
Vr0 = dot(r0, v0)/R0;

for i = 1:length(dnu)
    % FINDING NEXT POSITION VECTOR RADIUS
    Rf = (h^2/mu)*( 1 + ((h^2)/(mu*R0) - 1)*cos(dnu(i)) - ((h*Vr0)/mu)*sin(dnu(i)) )^-1;


    % CALCULATING f AND g
    f = 1 - ((mu*Rf)/h^2)*(1 - cos(dnu(i)));
    g = ((Rf*R0)/h)*sin(dnu(i));
    fdot = (mu/h)*( (Vr0/h)*(1 - cos(dnu(i))) - sin(dnu(i))/R0 );
    gdot = 1 - ((mu*R0)/h^2)*(1 - cos(dnu(i)));

    % CALCULATING VECTORS
    rf(i, :) = f*r0 + g*v0;
    vf(i, :) = fdot*r0 + gdot*v0;
end

% HANDLING SINGLE CASE
if length(dnu) == 1
    rf = rf.';
    vf = vf.';
end
end