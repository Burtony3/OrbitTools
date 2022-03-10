function [r, coe] = bodyPositionLF(id, jd, distU, angU)

if nargin < 3; distU = 'km'; end
if nargin < 4; angU = 'deg'; end

if strcmpi(distU, 'AU')
    DU = 1;
elseif strcmpi(distU, 'km')
    DU = 1 / 6.6845871226706E-9;
else
    DU = 1;
end

if strcmpi(angU, 'deg')
    AngU = 1;
else
    AngU = pi/180;
end

%% DATA TABLE
rates = struct;

% EARTH
id_ = 3;
rates(id_).a0  =     1.00000261;
rates(id_).da  =     5.62e-6;
rates(id_).e0  =     0.01671123;
rates(id_).de  =    -4.392e-5;
rates(id_).i0  =    -1.531e-5;
rates(id_).di  =    -0.01294668;
rates(id_).L0  =   100.46457166;
rates(id_).dL  = 35999.37244981;
rates(id_).om0 =   102.93768193;
rates(id_).dom =     0.32327364;
rates(id_).Om0 =     0;
rates(id_).dOm =     0;

% MARS
id_ = 4;
rates(id_).a0  =     1.53271034;
rates(id_).da  =     1.847e-5;
rates(id_).e0  =     0.09339410;
rates(id_).de  =     7.882e-5;
rates(id_).i0  =     1.84969142;
rates(id_).di  =    -0.00813131;
rates(id_).L0  =    -4.55343205;
rates(id_).dL  = 19140.30268499;
rates(id_).om0 =   -23.94362959;
rates(id_).dom =     0.44441088;
rates(id_).Om0 =    49.55953891;
rates(id_).dOm =    -0.29257343;

%% PULLING TABLE DATA
% FINDING CENTURIES PAST J2000
T = (jd - 2451545.0)/36525;

% GETTING ELEMENTS
a  = rates(id).a0 + rates(id).da*T;
e  = rates(id).e0 + rates(id).de*T;
i  = rates(id).i0 + rates(id).di*T;
L  = rates(id).L0 + rates(id).dL*T;
ombar = rates(id).om0 + rates(id).dom*T;
Om = rates(id).Om0 + rates(id).dOm*T;

% CONVERTING
om = ombar - Om;
M = L - ombar;
E = rad2deg(elm.M2E(e, deg2rad(M)));
coe = [a*DU e i*AngU Om*AngU om*AngU M*AngU];

% CONVERTING TO ECLIPTIC POSITION COORDINATES
xp = a*(cosd(E) - e);
yp = a*sqrt(1 - e^2)*sind(E);

r = [(cosd(om)*cosd(Om) - sind(om)*sind(Om)*cosd(i))*xp + (-sind(om)*cosd(Om) - cosd(om)*sind(Om)*cosd(i))*yp; 
     (cosd(om)*sind(Om) + sind(om)*cosd(Om)*cosd(i))*xp + (-sind(om)*sind(Om) + cosd(om)*cosd(Om)*cosd(i))*yp; 
                                  (sind(om)*sind(i))*xp +                               (cosd(om)*sind(i))*yp]*DU;
end