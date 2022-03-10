function [rSun, angSun] = solarPosition(date)

n = juliandate(date) - 2451545;
L = 280.46 + 0.9856474*n;
g = 357.528 + 0.9856003*n;

lambda = L + 1.915*sind(g) + 0.020*sind(2*g);
% beta = 0;
epsilon = 23.439 - 0.0000004*n;

rSun = [cosd(lambda); cosd(epsilon)*sind(lambda); sind(epsilon)*sind(lambda)];
angSun = atan2d(rSun(2), rSun(1));

end