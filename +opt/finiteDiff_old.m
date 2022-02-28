function g = finiteDiff_old(f, x, opts)
%FINITEDIFF finds the derivative of function f by taking (n-1) steps around
%   the center point x. For n=3, the points used would be f(x-h), f(x), and
%   f(x+h).
%
%   <strong>Function Calls:</strong>
%       g = FINITEDIFF(f, x)
%       g = FINITEDIFF(f, x, h=1e-8, n=5, usePar=true, ddx=2)
%
%   <strong>Required Inputs:</strong>
%       - f::<strong>Function Handle</strong>: Function of input vector x
%       - x::<strong>Vector{Double}</strong>: Point which derivative is
%           taken about
%
%   <strong>Name-Value Pairs:</strong>
%       - h::<strong>Double</strong>: Step length for each step.
%           Values between 1e-4 -- 1e-30 are valid
%           [Default] 1e-4
%       - n::<strong>Int</strong>: Number of points in stencil
%           [Default] 3
%       - ddx::<strong>Int</strong>: Which derivative to evaluate
%           [Default] 1
%       - usePar::<strong>Boolean</strong>: Parallelize derivative evaluation
%           [Default] false
%
%   <strong>Outputs:</strong>
%       - g::<strong>Vector/Matrix{Double}</strong>: Gradient of function f
%           at point x. Size nx x nf

% ARGUMENT VALIDATION
arguments
    f {mustBeSingleInputFunc}
    x (:, 1) {mustBeVector, mustBeNumeric}
    opts.h (1, 1) {mustBeNumeric, mustBePositive} = 1e-4
    opts.usePar (1, 1) {mustBeNumericOrLogical} = false
    opts.n (1, 1) {mustBeNumeric, mustBeScalarOrEmpty, mustBeInteger} = 3
    opts.ddx (1, 1) {mustBeNumeric, mustBeScalarOrEmpty, mustBeInteger} = 1
end
h = opts.h;
if opts.usePar; workers = Inf; else; workers = 0; end

% CREATING STENCILS
switch opts.n
    case 3
        A = [-0.5 0 0.5].';
        steps = -1:1;
    case 5
        A = [1/12, -2/3, 0, 2/3, -1/12].';
        steps = -2:2;
    case 7
        A = [-1/60, 9/60, -45/60, 0, 45/60, -9/60, 1/60].';
        steps = -3:3;
    otherwise
        A = [-0.5 0 0.5].';
        steps = -1:1;
end

% PREALLOCATING g
fval0 = f(x);
fsz = size(fval0);
if all(fsz == 1)
    g = zeros(length(x), 1);
elseif any(fsz == 1)
    g = zeros(length(x), fsz(fsz ~= 1));
else
    g = zeros(length(x), fsz(1), fsz(2));
end

% LOOPING THROUGH VALUES OF X
xinit = x;
for i = 1:length(x)
% parfor (i = 1:length(x), workers)
    % PREALLOCATING
    fval = [zeros(fsz(1), (opts.n-1)/2), fval0, zeros(fsz(1), (opts.n-1)/2)];

    % CREATING STEP
    x = xinit;
    H = zeros(length(x), 1);
    H(i) = h;

    % EVALUATING FUNCTION
    for j = 1:opts.n
        if A(j) == 0; continue; end
        fval(:, j) = f(x + steps(j)*H);
    end

    % CREATING g VARIABLE
%     if i == 1; g = zeros(length(x), length(fval(1).val)); end

    % SAVING TO VARIABLE
    g(i, :) = sum(A.*fval.', 1)./h;

end % END FOR

end % END FUNCTION

% VALIDATION FUNCTION
function mustBeSingleInputFunc(a)
    if ~isa(a, 'function_handle') && nargin(a) == 1
        eidType = 'mustBeSingleInputFunc:notBeSingleInputFunc';
        msgType = 'Input must be a single input, anonymous function.';
        throwAsCaller(MException(eidType,msgType))
    end
end

