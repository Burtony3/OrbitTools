function dX = CR3BP(~, X, mu)
% CR3BP is the equations of motion to propagate the state (X) over time 
%   interval (t), under the Circular Restrict Three Body Problem (rotating).
% 
% <strong>Function Calls:</strong>
%   - dX = THIRDBODY(t, X, muP, rP)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state. [Unused]
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time
%   - muP::<strong>Double</strong>: Stnd. grav. const. of third body
%   - rP::<strong>Vector{Double}</strong>: Radius from coordinate origin to perturbing 
%       body
% 
% <strong>Outputs:</strong>
%   - dX::<strong>Vector{Double}</strong>: Derivative of input vector at 
%       propagation time

% PULLING RADIUS VECTOR FROM STATE
error('Incomplete')
r = X(1:3);

% FINDING S/C -> PERTURBING BODY RADIUS
r_rP = r - rP;

% CALCULATING PERTURBATION
dX = [zeros(3, 1);
      -muP*( r_rP/norm(r_rP)^3 + rP/norm(rP)^3 )];

end