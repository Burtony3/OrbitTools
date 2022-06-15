function EOP = readEOP(utc, debug)
% READEOP collects and organizes coefficients used for ITRF->GCRF
%   transformation. Using C04 with Celestial Pole offsets (dPsi,dEps) 
%   referred to IAU 1980 precession-nutation model. 
% 
% <strong>Function Calls:</strong>
%   - EOP = READEOP(utc)
% 
% <strong>Required Inputs:</strong>
%   - utc::<strong>Vector{Datetime}</strong>: Array of dates to check. If there are 
%       two elements, READEOP will output the per date data in EOP.data and
%       linearly interopolated functions in EOP.func. Otherwise READEOP
%       will output just the data at the corresponding dates. 
% 
% <strong>Outputs:</strong>
%   - EOP::<strong>Struct</strong>: Output struct
%       - EOP.mjd::<strong>Vector{Double}</strong>: Modified Julian Date
%       - EOP.deltaUT1::<strong>Vector{Double}</strong>: Difference between UTC and UT1
%       - EOP.pmx::<strong>Vector{Double}</strong>: Polar motion in x-direction
%       - EOP.pmy::<strong>Vector{Double}</strong>: Polar motion in y-direction
%       - EOP.dPsi::<strong>Vector{Double}</strong>: Delta psi for nutation
%       - EOP.dEps::<strong>Vector{Double}</strong>: Delta epsilon for nutation

if nargin < 2; debug = false; end

% SETUP
utc2mjd = @(utc) floor(juliandate(utc, 'modifiedjuliandate'));

% GATHERING DATA
% https://hpiers.obspm.fr/eop-pc/index.php?index=C04&lang=en
site = urlread("https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.62-now");
site = strsplit(site, '\n').';
site = site(13:end-1);
tmp = str2num(char(site{1}));
data = zeros(length(site), length(tmp));
for i = 1:length(site)
    data(i, :) = str2num(char(site{i}));
end

% FINDING DATE BOUNDS
mjdbounds = [data(1, 4) data(end, 4)];
maxdatestring = sprintf('Limits are [%s, %s]', ...
    datestr(datetime(data(1, 1:3))), datestr(datetime(data(end, 1:3))) );

if debug; disp(site); end
usecase = length(utc);

%% COLLECTING VALUES
if usecase == 1     % Single-value output
    % CONVERTING JDATE
    mjd = utc2mjd(utc);
    if mjd < mjdbounds(1) || mjd > mjdbounds(2)
        error('Date pulled is out of bounds')
    end
    
    % PULLING FROM DATA
    out = collectvalues(data, mjd, debug);
    EOP = struct('mjd', mjd, 'deltaUT1', out.deltaUT1, ...
        'pmx', out.polarMotion(1), 'pmy', out.polarMotion(2), ...
        'dPsi', out.dPsi, 'dEps', out.dPsi);
elseif usecase == 2 % Range-output w/ curve-fits
    % CONVERTING TO MJD
    mjdrange = utc2mjd(utc);
    if any(mjdrange < mjdbounds(1) | mjdrange > mjdbounds(2))
        error('Date pulled is out of bounds\n%s', maxdatestring)
    end
    mjdrange = mjdrange(1):mjdrange(2);
    
    % PULLING DATA
    len = length(mjdrange);
    [EOP.data.mjd, EOP.data.deltaUT1, ...
     EOP.data.pmx, EOP.data.pmy, ...
     EOP.data.dPsi, EOP.data.dEps] = deal(zeros(len, 1));
    for i = 1:len
        out = collectvalues(data, mjdrange(i), debug);
        EOP.data.mjd(i) = mjdrange(i);
        EOP.data.deltaUT1(i) = out.deltaUT1;
        EOP.data.pmx(i) = out.polarMotion(1);
        EOP.data.pmy(i) = out.polarMotion(2);
        EOP.data.dPsi(i) = out.dPsi;
        EOP.data.dEps(i) = out.dEps;
    end
    
    % CREATING POLYNOMIALS
    fd = fields(EOP.data);
    for i = 1:length(fd)
        f = fit(EOP.data.mjd, EOP.data.(fd{i}), 'linearinterp');
        EOP.func.(fd{i}) = @(mjd) f(mjd);
    end
else                % Multi-value output
    % CONVERTING TO MJD
    mjdrange = utc2mjd(utc);
    if any(mjdrange < mjdbounds(1) | mjdrange > mjdbounds(2))
        error('Date pulled is out of bounds\n%s', maxdatestring)
    end
    
    % COLLECTING VALUES
    len = length(mjdrange);
    [EOP.mjd, EOP.deltaUT1, ...
     EOP.pmx, EOP.pmy, ...
     EOP.dPsi, EOP.dEps] = deal(zeros(len, 1));
    for i = 1:len
        out = collectvalues(data, mjdrange(i), debug);
        EOP.mjd(i) = mjdrange(i);
        EOP.deltaUT1(i) = out.deltaUT1;
        EOP.pmx(i) = out.polarMotion(1);
        EOP.pmy(i) = out.polarMotion(2);
        EOP.dPsi(i) = out.dPsi;
        EOP.dEps(i) = out.dEps;
    end
end

%% PULLS INFORMATION FROM MAIN DATA VARIABLE
function out = collectvalues(data, mjd, debug)
    % PULLING FROM CORRECT ROW
    idx = find(data(:, 4) == mjd, 1, 'first');
    data = data(idx, :);
    
    if debug; disp(data); end

    % COLLECTING INTO STRUCT
    out = struct;
    out.deltaUT1 = data(7);
    out.polarMotion = deg2rad(data(5:6)/3600);
    out.dPsi = deg2rad(data(9)/3600);
    out.dEps = deg2rad(data(10)/3600);
end
end