function r = rSunECI(utc, units, jd0)
% RSUNECI approximates the position of the Sun in the Earth-Centered
%   Inertial frame for a given date.
%
% From Vallado Algorithm 29
% 
% <strong>Function Calls:</strong>
%   - r = rSunECI(utc)
%   - r = rSunECI(utc, "km")
%   - r = rSunECI(seconds, "km", jd0)
% 
% <strong>Required Inputs:</strong>
%   - utc::<strong>Datetime/Double</strong>: Can either be a datetime input
%       <strong>or</strong> double of the Julian Date for the same datetime
%   - units::<strong>String</strong>: Output units, either "AU" or "km"
%   - jd0::<strong>Double</strong>: Initial Julian Date, where utc is the
%       measured seconds since this time
% 
% <strong>Outputs:</strong>
%   - r::<strong>Vector{Double}</strong>: Solar radius vector for given time

% ARGUMENT VALIDATION
if nargin < 2; units = 'AU'; end
if isdatetime(utc)
    jd = juliandate(utc);
else
    if nargin == 3
        seconds_plus_jd0 = utc;
        jd = jd0 + seconds_plus_jd0/86400;
    else
        jd = utc;
    end
end

% FINDING CENTURIES PAST J2000
jd2000 = juliandate(datetime('2000-01-01 12:00:00'));
T = (jd - jd2000) / 36525;

% MEAN SOLAR LONGITUDE AND ANOMALY
lam_M = 280.46 + 36000.771*T;
M = 357.5291092 + 35999.05034*T; 

% LONGITUDE OF ECLIPTIC
lam = lam_M + 1.914666471*sind(M) + 0.019994643*sind(2*M);

% OBLIQUITY OF ECLIPTIC
eps = 23.439291 - 0.0130042*T;

% SOLAR RADIUS SCALAR
R = 1.000140612 - 0.016708617*cosd(M) - 0.000139589*cosd(2*M);

% CONVERTING TO VECTOR
r = R*[cosd(lam); 
       cosd(eps)*sind(lam);
       sind(eps)*sind(lam)];

% OUTPUTTING
if strncmpi(units, 'km', 1)
    AU2km = 149597870.7;
    r = r*AU2km;
end

end