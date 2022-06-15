function r = rMoonECI(utc, jd0)
% RSUNECI approximates the position of the Sun in the Earth-Centered
%   Inertial frame for a given date.
%
% From Vallado Algorithm 31
% 
% <strong>Function Calls:</strong>
%   - r = rSunECI(utc)
%   - r = rSunECI(seconds, jd0)
% 
% <strong>Required Inputs:</strong>
%   - utc::<strong>Datetime/Double</strong>: Can either be a datetime input
%       <strong>or</strong> double of the julian date for the same datetime
%   - jd0::<strong>Double</strong>: Initial Julian Date, where utc is the
%       measured seconds since this time
% 
% <strong>Outputs:</strong>
%   - r::<strong>Vector{Double}</strong>: Solar radius vector for given time

% ARGUMENT VALIDATION
if isdatetime(utc)
    jd = juliandate(utc);
else
    if nargin == 2
        seconds_plus_jd0 = utc;
        jd = jd0 + seconds_plus_jd0/86400;
    else
        jd = utc;
    end
end

% FINDING CENTURIES PAST J2000
jd2000 = juliandate(datetime('2000-01-01 12:00:00'));
T = (jd - jd2000) / 36525;

% FINDING LONGITUDE w.r.t ECLIPTIC
lam = 218.32 + 481267.8813*T + 6.29*sind(134.9 + 477198.85*T) ...
    - 1.27*sind(259.2 - 413335.38*T) + 0.66*sind(235.7 + 890534.23*T) ...
    + 0.21*sind(269.9 + 954397.70*T) - 0.19*sind(357.5 + 35999.05*T) ...
    - 0.11*sind(186.6 + 966404.05*T);

% FINDING LATITUDE w.r.t. ECLIPTIC
phi = 5.13*sind(93.3 + 483202.03*T) + 0.28*sind(228.2 + 960400.87*T) ...
    - 0.28*sind(318.3 + 6003.18*T) - 0.17*sind(217.6 - 407332.20*T);

% FINDING PARALLAX
p = 0.9508 + 0.0518*cosd(134.9 + 477198.85*T) ...
    + 0.0095*cosd(259.2 - 413335.38*T) + 0.0078*cosd(235.7 + 890534.23*T) ...
    + 0.0028*cosd(269.9 + 954397.70*T);

% ANGLE TO ECLIPTIC
eps = 23.439291 - 0.0130042*T - 1.64e-7*T^2 + 5.04e-7*T^3;

% FINDING DISTANCE USING PARALLAX
basis = 6378.137;                   % Radius of Earth
R = basis/sind(p);

% FINDING VECTOR
r = R*[cosd(phi)*cosd(lam);
       cosd(eps)*cosd(phi)*sind(lam) - sind(eps)*sind(phi);
       sind(eps)*cosd(phi)*sind(lam) + cosd(eps)*sind(phi)];

end