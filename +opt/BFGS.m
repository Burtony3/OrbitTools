function out = BFGS(x0, f, g, opts)
%BFGS navigates to the bottom of function, f, given the gradient, g,
%   starting from initial point, x0. Uses Broyden-Fletcher-Goldfarb-Shanno
%   algorithm to approximate hessian of f for more accurate gradient
%   descent direction choices.
%
%   <strong>Function Calls:</strong>
%       bfgs = BFGS(x0, f, g)
%       bfgs = BFGS(x0, f, g, opts)
%
%   <strong>Required Inputs:</strong>
%       - x0::<strong>Vector{Double}</strong>: Starting point for
%           optimization
%       - f::<strong>Function Handle</strong>: Function of input vector x
%       - g::<strong>Function Handle</strong>: Gradient of f with input vector x
%
%   <strong>Name-Value Pairs:</strong>
%       - Tol::<strong>Double</strong>: Optimization minimum tolerance
%           [Default] 1e-8
%       - IterBreak::<strong>Int</strong>: Maximum number of iterations
%           [Default] 3
%       - InitialStride::<strong>Double</strong>: Initial stride length for
%           line search.
%           [Default] 0.1
%       - StrideNum::<strong>Int</strong>: Number of strides before
%           increasing step size.
%           [Default] 4
%       - StrideMod::<strong>Double</strong>: Stride length multiplier
%           [Default] 2
%       - IterPlot::<strong>Function Handle</strong>: Function to create
%           plot each iteration with input of x
%           [Default] 0x0 Empty
%
%   <strong>Outputs:</strong>
%       - out::<strong>Struct</strong>: Collection of output variables
%           <strong>Fields:</strong>
%           - x::<strong>Vector{Double}</strong>: Minimum value state
%           - fCalls::<strong>Vector{Double}</strong>: Number of
%               Function and Gradient calls
%           - xhist::<strong>Matrix{Double}</strong>: Value of x at
%               each iteration
%           - iters::<strong>Int</strong>: Number of minor iterations

%% ARGUMENT VALIDATION
arguments
    x0 (:, 1) {mustBeVector, mustBeNumeric}
    f {mustBeSingleInputFunc}
    g {mustBeSingleInputFunc}
    opts.Tol (1, 1) {mustBeNumeric, mustBePositive} = 1e-8
    opts.IterBreak (1, 1) {mustBeNumeric, mustBePositive, mustBeInteger} = 1e3
    opts.InitialStride (1, 1) {mustBeNumeric, mustBePositive} = 0.1
    opts.StrideNum (1, 1) {mustBeNumeric, mustBePositive, mustBeInteger} = 4
    opts.StrideMod (1, 1) {mustBeNumeric, mustBePositive} = 2
    opts.IterPlot {mustBeSingleInputFunc}
    opts.UsingFD (1, 1) {mustBeNumericOrLogical} = false
end
tol = opts.Tol;
iterBreak = opts.IterBreak;
t0 = opts.InitialStride;
steps = opts.StrideNum;
stepsMod = opts.StrideMod;
if isfield(opts, 'IterPlot')
    IterPlot = opts.IterPlot;
    plotIters = true;
else
    plotIters = false;
end

%% STARTUP
dx = Inf;
i = 0;
dir0 = [];
fCalls = [0 0];
x = x0;
xhist = x0(:);
dxhist = [];

while dx > tol
    i = i+1;
    % FINDING DIRECTION
    dir = BFGS_Update(x, g, dir0, opts.UsingFD); % Broyden-Fletcher-Goldfarb-Shanno
    s = dir.s;
    fCalls = fCalls + dir.fCalls;

    % LINE SEARCHING
    lin = lineSearch(x, s, f, t0, steps, stepsMod);
%     fCalls = fCalls + lin.fCalls;

    % MINIMIZING
    mini = goldenRatioMin([lin.pts(:).x], f);
%     fCalls = fCalls + mini.fCalls;

    % ALIGNING NEW MINIMUM
    dx = norm(x(2:end) - mini.xmin(2:end));
    dxhist(end+1) = dx;
    x = mini.xmin;
    xhist = [xhist x(:)];

    % UPDATING FOR NEXT STEP
    dir0 = dir;
    if i == iterBreak; break; end

    % PLOTTING
    if plotIters; IterPlot(x); end
end % END WHILE

% OUTPUTTING
out.x = x;
out.fCalls = fCalls;
out.xhist = xhist;
out.iters = i;
out.dx = dxhist;

%% SUB-ROUTINE
function out = BFGS_Update(x, g, prev, FDBool)
    if isempty(prev)
        gk1 = g(x);
        Qk1 = eye(length(x));
    else
        % EVALUTING NEW TERMS
        gk1 = g(x);
        Qk = prev.Q;

        % CREATING CONSTANTS
        p = x - prev.x;
        y = gk1 - prev.g;
        sig = p.' * y;
        tau = y.' * Qk * y;
        A = Qk * y * p.';

        % COMPUTING UPDATE
        dQ = ((sig + tau)/sig^2)*(p * p.') - (1/sig)*(A + A.');
        Qk1 = Qk + dQ; % Creating approximated Hessian
    end % END IF

    out.s = -Qk1*gk1; out.s = out.s/norm(out.s);
    out.x = x;
    if FDBool
        out.fCalls = [length(x)*3 1];
    else
        out.fCalls = [0 1];
    end
    out.g = gk1;
    out.Q = Qk1;
end % END SUB-ROUTINE
end % END FUNCTION

%% VALIDATION FUNCTION
function mustBeSingleInputFunc(a)
    if ~isa(a, 'function_handle') && nargin(a) == 1
        eidType = 'mustBeSingleInputFunc:notBeSingleInputFunc';
        msgType = 'Input must be a single input, anonymous function.';
        throwAsCaller(MException(eidType,msgType))
    end
end
