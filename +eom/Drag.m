function dX = Drag(~, X, CD, A, m, rho0, r0, H, thetadot)
% DRAG is the equations of motion to propagate the state (X) over time 
%   interval (t), under the inertially-fixed, cannonball drag model of
%   spacecraft drag.
% 
% <strong>Function Calls:</strong>
%   - dX = DRAG(t, X, CD, A, m, rho0, r0, H, thetadot)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state. [Unused]
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time  (km & km/s)
%   - CD::<strong>Double</strong>: Coefficient of drag for spacecraft
%   - A::<strong>Double</strong>: Cross-sectional area normal to velocity  (m^2)
%   - m::<strong>Double</strong>: Mass of spacecraft                       (kg)
%   - rho0::<strong>Double</strong>: Atmospheric density at height of r0   (kg/m^3)
%   - r0::<strong>Double</strong>: See above                               (m)
%   - H::<strong>Double</strong>: Scale height for atmosphere model        (m)
%   - thetadot::<strong>Double</strong>: Rotational rate of planet         (1/sec)
% 
% <strong>Outputs:</strong>
%   - dX::<strong>Vector{Double}</strong>: Derivative of input vector at   (km/s^2)
%       propagation time

% CONSTANT
km2m = 1000;

% VELOCITY W.R.T. ATMOSPHERE
vvec = [X(4) + thetadot*X(2);
        X(5) - thetadot*X(1);
        X(6)];
v = sqrt(sum(vvec.^2));

% CALCULATING DENSITY
rho = rho0*exp(-(km2m*norm(X(1:3)) - r0)/H); % 1000 converts from km -> m

% CALCULATING ACCELERATION
dX = [zeros(3, 1); 
      -km2m*0.5*CD*(A/m)*rho*(v*vvec)]; % 1000 accounts for extra m from density
end