function state = meanState(utc, id, opts)


%% ARGUMENT VALIDATION
arguments
    utc datetime
    id (1, 1) {mustBeInteger}
    opts.ReturnElements = false
    opts.ReturnTableElements = false;
    opts.Frame {mustBeMember(opts.Frame, ["ECLIPJ2000", "J2000", "HELIO"])} = "ECLIPJ2000" 
end

%
lenout = length(utc);
jd2000 = juliandate(datetime('2000-01-01'));
if year(utc) > 2050
    error('Year must be less than 2050')
end

%% TABLES
%           a                e                i                 L            long.peri.       long.node.
%          AU               rad              deg               deg              deg              deg
elmt = [0.38709927,      0.20563593,      7.00497902,      252.25032350,     77.45779628,     48.33076593;   %mercury
        0.72333566,      0.00677672,      3.39467605,      181.97909950,    131.60246718,     76.67984255;   %venus
        1.00000261,      0.01671123,     -0.00001531,      100.46457166,    102.93768193,      0.0;          %earth moon barycenter
        1.52371034,      0.09339410,      1.84969142,       -4.55343205,    -23.94362959,     49.55953891;   %mars
        5.20288700,      0.04838624,      1.30439695,       34.39644051,     14.72847983,    100.47390909;   %jupiter
        9.53667594,      0.05386179,      2.48599187,       49.95424423,     92.59887831,    113.66242448;   %saturn
        19.18916464,     0.04725744,      0.77263783,      313.23810451,    170.95427630,     74.01692503;   %uranus
        30.06992276,     0.00859048,      1.77004347,      -55.12002969,     44.96476227,    131.78422574;   %neptune
        39.48211675,     0.24882730,     17.14001206,      238.92903833,    224.06891629,    110.30393684 ]; %pluto

%          AU/Cy           rad/Cy           deg/Cy           deg/Cy              deg/Cy           deg/Cy
rate = [ 0.00000037,      0.00001906,     -0.00594749,   149472.67411175,      0.16047689,     -0.12534081;  %mercury
         0.00000390,     -0.00004107,     -0.00078890,    58517.81538729,      0.00268329,     -0.27769418;  %venus
         0.00000562,     -0.00004392,     -0.01294668,    35999.37244981,      0.32327364,      0.0;         %earth moon barycenter
         0.00001847,      0.00007882,     -0.00813131,    19140.30268499,      0.44441088,     -0.29257343;  %mars
        -0.00011607,     -0.00013253,     -0.00183714,     3034.74612775,      0.21252668,      0.20469106;  %jupiter
        -0.00125060,     -0.00050991,      0.00193609,     1222.49362201,     -0.41897216,     -0.28867794;  %saturn
        -0.00196176,     -0.00004397,     -0.00242939,      428.48202785,      0.40805281,      0.04240589;  %uranus
         0.00026291,      0.00005105,      0.00035372,      218.45945325,     -0.32241464,     -0.00508664;  %neptune
        -0.00031596,      0.00005170,      0.00004818,      145.20780515,     -0.04062942,     -0.01183482]; %pluto
    
%% COMPUTING ELEMENTS
coe = zeros(lenout, 6);
for i = 1:lenout
    % CENTURIES PAST J2000
    jd = juliandate(utc(i));
    T = (jd - jd2000)/36525;
    
    % COMPUTING ELEMENTS
    for j = 1:6
        coe(i, j) = elmt(id, j) + rate(id, j)*T;
    end
end
if opts.ReturnTableElements
    state = coe;
    if lenout == 1; state = state.'; end
    return
end

% SWAPPING LAST ELEMENTS
om = coe(:, 5) - coe(:, 6);
M  = coe(:, 4) - coe(:, 5);
coe(:, 3) = mod(coe(:, 3), 360);
coe(:, 4) = mod(coe(:, 6), 360);
coe(:, 5) = mod(om, 360);
coe(:, 6) = mod(M, 360);

% BOUNDING ANGLES
% for i = 3:6
%     coe(:, i) = mod(coe(:, i), 360);
% end

% ENDING IF OUTPUTTING ELEMENTS
if opts.ReturnElements
    state = coe;
    if lenout == 1; state = state.'; end
    return
end

%% CONFIGURING OUTPUTS
% CALCULATING ECCENTRIC ANOMALY
E = zeros(lenout, 1);
for i = 1:lenout
    E(i) = rad2deg( elm.M2E( coe(i, 2), deg2rad(coe(i, 6)) ) );
end

% FINDING HELIOCENTRIC COORDINATES
r = [coe(:, 1).*( cosd(E) - coe(:, 2) ), coe(:, 1).*sqrt( 1 - coe(:, 2).^2 ).*sind(E), zeros(lenout, 1)];
if strcmpi(opts.Frame, "HELIO")
    state = r;
    if lenout == 1; state = state.'; end
    return
end

% FINDING ECLIPTIC COORDINATES
r_ecl = zeros(lenout, 3);
for i = 1:lenout
    r_ecl(i, :) = frame.rot('313', [-coe(i, 4) -coe(i, 3) -coe(i, 5)], 'deg')*r(i, :).';
end
if strcmpi(opts.Frame, "ECLIPJ2000")
    state = r_ecl;
    if lenout == 1; state = state.'; end
    return
end

% FINDING EQUITORIAL COORDINATES
r_j2000 = zeros(lenout, 3);
inc_j2000 = 23.43928;
R = frame.rot('1', -inc_j2000, 'deg');
for i = 1:lenout
    r_j2000(i, :) = R*r_ecl(i, :).';
end

state = r_j2000;
if lenout == 1; state = state.'; end
end