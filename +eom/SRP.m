function dX = SRP(~, X, CR, A, m, rS)
% SRP is the equations of motion to propagate the state (X) over time 
%   interval (t), under the inertially-fixed, sollar radiation pressure.
%
% From Vallado Section 8.6.4 pg. 578 (606)
% 
% <strong>Function Calls:</strong>
%   - dX = SRP(t, X, CR, A, m, rS)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state. [Unused]
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time  (km & km/s)
%   - CR::<strong>Double</strong>: Coefficient of reflectivity (0 - 2)
%   - A::<strong>Double</strong>: Affected area                            (m^2)
%   - m::<strong>Double</strong>: Spacecraft mass                          (kg)
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

% SCALING PRESSURE BY DISTANCE FROM REFERENCE
ratio = Rsc2S / AU;
scale = 1/ratio^2;

% TURNING RADIUS VECTOR TO UNIT
rsc2S = rsc2S / Rsc2S;

% ACCOUNTING FOR UNITS
m2km = 1/1000;          % Accounts for extra m in Psrp

% CONSTRUCTING ACCELERATION
dX = [zeros(3, 1);
      -m2km*scale*Psrp*CR*A/m * rsc2S]; 

end