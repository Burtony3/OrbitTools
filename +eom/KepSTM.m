function dX = KepSTM(~, X, mu)
% KEP is the equations of motion to propagate the state (X) over time 
%   interval (t), under the inertially-fixed Keplarian graviational force 
%   model. State includes 6x6 state transition matrix.
% 
% <strong>Function Calls:</strong>
%   - dX = KEPSTM(t, X, mu)
% 
% <strong>Required Inputs:</strong>
%   - t::<strong>Double</strong>: Current time of propagation related to 
%       the state
%   - X::<strong>Vector{Double}</strong>: State vector at propgation time
%   - mu::<strong>Double</strong>: Standard Gravitational Parameter
% 
% <strong>Outputs:</strong>
%   - dX::<strong>Vector{Double}</strong>: Derivative of input vector at 
%       propagation time

    % PULLING VALUES FROM STATE
    r = sqrt(sum(X(1:3).^2));
    R = X(1:3);
    phi = reshape(X(7:42), 6, 6);
    
    % CALCULATING STM INFORMATION
    G = mu/r^5 * ( 3*(R*R.') - (r^2 * eye(3)) );  % From BATES
    A = [zeros(3) eye(3); G zeros(3)];

    % DEFINING DERIVATIVE
    dX = [X(4:6);
        -mu/r^3*R;
        reshape(A*phi, 36, 1)];
end