function g = cplxDiff(f, x, h)
%CPLXDERIVATIVE finds the derivative of function f by taking a step in the
%   complex plane for each index of x with a step length of h.
%
%   <strong>Function Calls:</strong>
%       g = CPLXDERIVATIVE(f, x, h)
%
%   <strong>Inputs:</strong>
%       - f::<strong>Function Handle</strong>: Function of input vector x
%       - x::<strong>Vector{Double}</strong>: Point which derivative is
%           taken about
%       - h::<strong>Double</strong>: Step length for each step.
%           Values between 1e-4 -- 1e-30 are valid
%
%   <strong>Outputs:</strong>
%       - g::<strong>Vector/Matrix{Double}</strong>: Gradient of function f
%           at point x

% LOOPING THROUGH VALUES OF X
xinit = x;
for i = 1:length(x)

    % CREATING STEP
    x = xinit;
    x(i) = x(i) + h*1i;

    % EVALUATING FUNCTION
    fval = f(x);

    % CREATING g VARIABLE
    if i == 1; g = zeros(length(x), length(fval)); end

    % SAVING TO VARIABLE
    g(i, :) = imag(fval)/h;

end % END FOR

end % END FUNCTION