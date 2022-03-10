function varargout = nondimension(rvec, vvec, mu)
% NONDIMM non-dimensionalizes the system given the input state. Converting
%   positions and velocities to Distance Units and Time Units, based upon a
%   circular orbit at the radius of ||rvec||. Outputted system parameters
%   will have a mu = 1, and a ||rvec|| = 1. 
% 
% <strong>Function Calls:</strong>
%   - [DU, TU] = NONDIMM(rvec, vvec, mu)
%   - [DU, TU, rvec, vvec] = NONDIMM(rvec, vvec, mu)
% 
% <strong>Required Inputs:</strong>
%   - rvec::<strong>Vector{Double}</strong>: Orbital Position Vector
%   - vvec::<strong>Vector{Double}</strong>: Orbital Velocity Vector
%   - mu::<strong>Double</strong>: Standard Gravitational Parameter
% 
% <strong>Outputs:</strong>
%   - DU::<strong>Double</strong>: Distance units conversion to inputted
%       units
%   - TU::<strong>Double</strong>: Time units conversion to inputted units
%   - rvec::<strong>Vector{Double}</strong>: Non-Dimensional Orbital Position Vector
%   - vvec::<strong>Vector{Double}</strong>: Non-Dimensional Orbital Velocity Vector

% NON-DIMENSIONAL DISTANCE UNITS
DU = norm(rvec);                % km / DU

% NON-DIMENSIONAL TIME UNITS
vc = sqrt(mu/DU);
TU = 1/(vc/DU);                 % seconds / TU

% CONVERTING UNITS
rvec = rvec / DU;
vvec = vvec / (DU / TU); 

% OUTPUTTING
if nargout == 2
    varargout = {DU, TU};
elseif nargout == 4
    varargout = {DU, TU, rvec, vvec};
else
    error('Number of outputs do not match function optionals')
end

end