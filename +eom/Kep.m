function dX = Kep(~, X, mu)
% KEP is the equations of motion to propagate the state (X) over time 
%   interval (t), under the inertially-fixed Keplarian graviational force 
%   model.
% 
% <strong>Function Calls:</strong>
%   - dX = KEP(t, X, mu)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time
%   - mu::<strong>Double</strong>: Standard Gravitational Parameter
% 
% </strong>Outputs:</strong>
%   - dX::<strong>Vector{Double}</strong>: Derivative of input vector at 
%       propagation time

    % PULLING VALUES FROM STATE
    r = sqrt(sum(X(1:3).^2));
    R = X(1:3);

    % DEFINING DERIVATIVE
    dX = [X(4:6);
        -mu/r^3*R];
end