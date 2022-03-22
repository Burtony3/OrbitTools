function [r, v, R] = inertial2perifocal(r, v, mu)

%% ARGUMENT VALIDATION
r = r(:);
v = v(:);

%% FINDING FOUNDATIONAL VECTORS
% h = cross(r, v);
% w = h/norm(h);
% 
% e = cross(v, h)/mu - r/norm(r);
% p = e/norm(e);
% 
% q = cross(w, p);

%% FINDING ROTATION MATRIX
% GETTING ORBITAL ELEMENTS
coe = elm.rv2coe(r, v, Mu = mu, Anomaly = 'true', Units = 'deg');
Omega = coe(4);
omega = coe(5);
i = coe(3);

% CALCULATING FRAME
R = frame.rot('313', [-Omega -i -omega], 'deg');

%% ROTATING AND REMOVING ELEMENTS
r = R.'*r; r = r(1:2);
v = R.'*v; v = v(1:2);

end