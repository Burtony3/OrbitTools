function [dOmega, domega, dM] = precessionJ2(a, e, i, mu, J2, ae)
% PRECESSIONJ2 generates rates of change for orbital elements under
%   spherical harmonics, given initial orbit parameters and system
%   dynamics.
% 
% <strong>Function Calls:</strong>
%   - [dΩ, dω, dM] = PRECESSIONJ2(a, e, i, μ, J2, ae)
% 
% <strong>Required Inputs:</strong>
%   - a::<strong>Double</strong>: Semi-Major Axis
%   - e::<strong>Double</strong>: Eccentricity
%   - i::<strong>Double</strong>: Inclination
%   - μ::<strong>Double</strong>: Standard Gravitational Parameter
%   - J2::<strong>Double</strong>: Oblateness Factor
%   - ae::<strong>Double</strong>: Reference Distance
% 
% <strong>Outputs:</strong>
%   - dΩ::<strong>Double</strong>: Rate of Change of LAAN
%   - dω::<strong>Double</strong>: Rate of Change of Argument of Periapsis
%   - dM::<strong>Double</strong>: Rate of Change of Mean Anomaly

% USEFUL CONSTANTS
n = sqrt(mu/a^3);
J0 = 1.5*n*J2*(ae/a)^2;

% CHANGE IN LAAN
dOmega = -J0/(1 - e^2)^2 * cosd(i)*180/pi;

% CHANGE IN ARGUMENT OF PERIAPSIS
domega = -J0/(1 - e^2)^2*(2.5*sind(i)^2 - 2)*180/pi;

% CHANGE IN MEAN ANOMALY
dM = n - J0/(1 - e^2)^(3/2) * (0.5 - 1.5*cosd(i)^2)*180/pi;
end