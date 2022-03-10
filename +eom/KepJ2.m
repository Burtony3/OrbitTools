function dX = KepJ2(~, X, mu, J2, ae)
% KEPJ2 is the equations of motion to propagate the state (X) over time 
%   interval (t), under the inertially-fixed, spherical harmonics 
%   graviational force model.
% 
% <strong>Function Calls:</strong>
%   - dX = KEPJ2(t, X, mu, J2, ae)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time
%   - mu::<strong>Double</strong>: Standard Gravitational Parameter
%   - J2::<strong>Double</strong>: Oblateness Factor
%   - ae::<strong>Double</strong>: Reference Distance
% 
% <strong>Outputs:</strong>
%   - dX::<strong>Vector{Double}</strong>: Derivative of input vector at 
%       propagation time

    % PULLING VALUES FROM STATE
    r = sqrt(sum(X(1:3).^2));
    R = X(1:3);

    % GENERATING CONSTANTS
    J0 = 1.5*mu*J2*ae^2;
    Jr = (J0 / r^5)*(1 - (5 / r^2)*dot(R, [0; 0; 1])^2);
    Jk = 2*(J0 / r^5)*dot(R, [0; 0; 1]);
    P = Jr*R + Jk*[0; 0; 1];

    % DEFINING DERIVATIVE
    dX = [X(4:6);
        -mu/r^3*R - P];
end