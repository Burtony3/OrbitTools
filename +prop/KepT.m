function X = KepT(rvec, vvec, tspan, mu)
%KEPT Summary of this function goes here
%   Detailed explanation goes here
%
%   X = KEPT(rvec, vvec, tspan, mu)

% ASSIGNING CONSTANTS
r = norm(rvec); v = norm(vvec);
% mu = env.mu;
eps = (0.5*v^2 - mu/r);
a = -0.5*mu/eps;              %
sig = dot(rvec,vvec)/sqrt(mu);           % Frequent Constant

if (isstring(tspan) || ischar(tspan)) && strcmpi(tspan, 'period')
    T = 2*pi*sqrt(a^3 / mu);
    tspan = linspace(0, T, 100);
end

[rvec_out, vvec_out] = deal(zeros(length(tspan), 3));

for i = 1:length(tspan)
    % FINDING CONIC TYPE
    if a > 0
        % DEFINING DIFFERENTIAL EQUATION AND FINDING ROOT
        dEi = tspan(i)*sqrt(mu/a^3) ;                                                 % Initial Guess for dE
        fdE = @(dE) (-dEi + dE + (sig/sqrt(a))*(1 - cos(dE)) - (1 - r/a)*sin(dE));    % Pykep kepler_equations.hpp line 49 
        dE = fsolve(fdE, dEi, optimset('Display', 'off'));

        % CREATING USEFUL CONSTANT
        rho = a + (r - a)*cos(dE) + sig*(sqrt(a))*sin(dE);                             % Battin Eqn. 4.42

        % FINDING LAGRANGIAN COEFFICIENTS
        F = 1 - (a / r) * (1 - cos(dE));                                      % Battin Eqn. 3.42
        G = a * sig / sqrt(mu) * (1 - cos(dE)) + r * sqrt(a / mu) * sin(dE);           % ""
        Ft = -sqrt(mu * a) / (rho * r) * sin(dE);                                % ""
        Gt = 1 - a / rho * (1 - cos(dE));                                       % ""

    else
        if tspan(i) > 0
            psi = 1; 
        else
            psi = -1;
        end
        disp('Hyperbolic')
        
%         dH = (rvec, ψ, rho, sig, a) -> (-psi - rho + (sig/√a)*(cosh(rho) - 1) + (1 - norm(rvec)/a)*sinh(rho))
    end

    % SOLVING FOR NEW rvec & vvec VECTORS
    rvec_out(i, :) = F*rvec + G*vvec;                  % Battin Eqn. 3.33
    vvec_out(i, :) = Ft*rvec + Gt*vvec;                % ""
end

X = [rvec_out vvec_out];

end

