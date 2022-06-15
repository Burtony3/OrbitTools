function dX = SRP_SD(~, X, Cd, Cs, A, m, rBody, rS)
% SRP_SD is the equations of motion to propagate the state (X) over time 
%   interval (t), under the inertially-fixed, sollar radiation pressure.
%
% From Vallado Section 8.6.4 pg. 582 (610)
% 
% <strong>Function Calls:</strong>
%   - dX = SRP_SD(t, X, Cd, Cs, A, m, rS)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state. [Unused]
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time  (km & km/s)
%   - Cd::<strong>Double</strong>: Diffusive reflection coefficient        (ND)
%   - Cs::<strong>Double</strong>: Specular reflection coefficient         (ND)
%   - A::<strong>Double</strong>: Affected area                            (m^2)
%   - m::<strong>Double</strong>: Spacecraft mass                          (kg)
%   - rBody::<strong>Vector{Double}</strong>: Unit vector normal to face   (ND)
%   - rS::<strong>Vector{Double}</strong>: Radius from coordinate origin   (km)
%       to Sun
% 
% <strong>Outputs:</strong>
%   - dX::<strong>Vector{Double}</strong>: Derivative of input vector at 
%       propagation time

% PRESSURE AT 1 AU 
AU = 149597871;
Psrp = 4.57e-6; % N/m^2 (For absorbative material, 2x for reflective)

% VECTOR FROM SATELLITE TO SUN
rsc2S = X(1:3) - rS;
Rsc2S = sqrt(sum(rsc2S.^2));
rsc2S = rsc2S / Rsc2S;          % Unit Vectorizing

% SCALING PRESSURE BY DISTANCE FROM REFERENCE
ratio = Rsc2S / AU;
scale = 1/ratio^2;
Psrp = scale*Psrp;

% GETTING UNIT VECTORS
n = rBody / norm(rBody);
s = -rsc2S;
phi = abs(frame.u2vAng(s, n, 'd', 'inside')); 

% CALCULATING SPECULAR FORCE
frs = -2*Psrp*Cs*A*cosd(phi)^2*n;

% CALCULATING DIFFUSIVE FORCE
frd = -Psrp*Cd*A*cosd(phi)*((2/3)*n + s);

% ACCOUNTING FOR UNITS
m2km = 1/1000;          % Accounts for extra m in Psrp

% CONSTRUCTING ACCELERATION
dX = [zeros(3, 1);
      m2km*(frs + frd)/m]; 


end