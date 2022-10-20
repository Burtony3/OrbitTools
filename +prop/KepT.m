function X = KepT(rvec, vvec, tspan, mu)
%KEPT Summary of this function goes here
%   Detailed explanation goes here
%
%   X = KEPT(rvec, vvec, tspan, mu)
%   tspan is in seconds

%% SETUP
% ASSIGNING CONSTANTS
r = norm(rvec); v = norm(vvec);
eps = (0.5*v^2 - mu/r);
a = -0.5*mu/eps;              %
sig = dot(rvec,vvec)/sqrt(mu);           % Frequent Constant

% DEFINING FUNCTIONS DEPENDING ON ORBIT
if a > 0
    % FINDING ECCENTRIC ANOMALY
    fxi = @(t) t*sqrt(mu/a^3);                                                % Initial Guess for dE
    f = @(x, t) (-sqrt(mu / a^3)*t + x + (sig/sqrt(a))*(1 - cos(x)) - (1 - r/a)*sin(x));   % Pykep kepler_equations.hpp line 49 

    % USEFUL CONSTANT
    frho = @(x) a + (r - a)*cos(x) + sig*(sqrt(a))*sin(x);                    % Battin Eqn. 4.42

    % F AND G
    fF  = @(x) 1 - (a / r) * (1 - cos(x));                                    % Battin Eqn. 3.42
    fG  = @(x) a * sig / sqrt(mu) * (1 - cos(x)) + r * sqrt(a / mu) * sin(x); % ""
    fdF = @(x, rho) -sqrt(mu * a) / (rho * r) * sin(x);                       % ""
    fdG = @(x, rho) 1 - a / rho * (1 - cos(x));                               % ""

else
    % FINDING ECCENTRIC ANOMALY
    fxi = @(t) sign(t);                                                       % Initial Guess for dH
    f = @(x, t) (-sqrt(-mu / a^3)*t - x + sig / sqrt(-a) * (cosh(x) - 1) + (1 - norm(rvec) / a) * sinh(x));

    % USEFUL CONSTANT
    frho = @(x) a + (r - a)*cosh(x) + sig*sqrt(-a)*sinh(x);                   % Battin Eqn. 4.42

    % F AND G
    fF  = @(x) 1 - a/r*(1 - cosh(x));
    fG  = @(x) a*sig/sqrt(mu)*(1 - cosh(x)) + r*sqrt(-a/mu)*sinh(x);
    fdF = @(x, rho) -sqrt(-mu*a)/(rho*r)*sinh(x);
    fdG = @(x, rho) 1 - a/rho*(1 - cosh(x));

end

% HANDLING SPECIAL CASE
if (isstring(tspan) || ischar(tspan)) && strcmpi(tspan, 'period')
    T = 2*pi*sqrt(a^3 / mu);
    tspan = linspace(0, T, 100);
    if a < 0; error('Orbit is hyperbolic: ''period'' is not allowed'); end
end

% PREALLOCATING
[rvec_out, vvec_out] = deal(zeros(length(tspan), 3));

%% PROPAGATING
for i = 1:length(tspan)
    % ROOT SOLVING FROM INITIAL CONDITIONS
    xi = fxi(tspan(i));
    xstar = fzero(@(x) f(x, tspan(i)), xi, optimset('Display', 'off'));

    % CREATING CONSTANT
    rho = frho(xstar);

    % FINDING F and G
    F = fF(xstar);
    G = fG(xstar);
    dF = fdF(xstar, rho);
    dG = fdG(xstar, rho);

    % SOLVING FOR NEW rvec & vvec VECTORS
    rvec_out(i, :) = F*rvec + G*vvec;                  % Battin Eqn. 3.33
    vvec_out(i, :) = dF*rvec + dG*vvec;                % ""
end

%% OUTPUTTING
X = [rvec_out vvec_out];

end



