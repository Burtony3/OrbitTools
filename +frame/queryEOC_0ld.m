function s = queryEOC(tspan, step)
% itrf2icrf queries Earth Orientation Center website for a given timespan
% and stepsize, returning Earth Orientation rotation matrix from
% Iternational Terrestrial Reference Frame to International Celestial
% Reference Frame.
%
% <strong>Warning</strong>: ~1 second run time, suggested to pull range of
% values and interpolate results on the fly.
%
% Transformation from the international terrestrial reference system (ITRF) 
% to the international celestial reference frame (ICRF): transformation 
% coordinate and associated quaternion. Compution of Celestial Pole 
% coordinates (X,Y), Earth rotation angle, Celestial Intermediate Origin 
% locator, Terrestrial Intermediate Origin locator are performed by FORTRAN 
% programs of the SOFA IAU library, consistent with the IAU 2000/2006 
% resolutions.
%
% EOC Link: https://hpiers.obspm.fr/eop-pc/index.php?index=rotation&lang=en
% 
% <strong>Function Calls:</strong>
%   - s = Ritrf2icrf(tspan, step)
% 
% <strong>Required Inputs:</strong>
%   - tspan::<strong>Vector{Datetime}</strong>: 2 element vector of dates
%   - step::<strong>Integer</strong>: Seconds between time querries in tspan
% 
% <strong>Optional Inputs:</strong>
%   - n/a
% 
% <strong>Outputs:</strong>
%   - s::<strong>Vector{Struct}</strong>: Struct containing the following fields:
%       - MJD::<strong>Double</strong>: Julian Date of query entry
%       - Date::<strong>Datetime</strong>: Datetime of query entry
%       - R::<strong>Matrix{Double}</strong>: 3x3 Rotation Matrix of entry
%       - Quat::<strong>Vector{Double}</strong>: Quaternion pointing vector
%           of Earth Orientation
% 
% Also see: 

% ARGUMENT VALIDATION
if nargin < 2; step = 3600; end
if ~isdatetime(tspan); error('Input 1 "tspan" must be a datetime vector'); end

% URL SETUP
url = 'https://hpiers.obspm.fr/eop-pc/products/matrice/rotation.php?';
params = ['&nut=1&ut1=1&pm=1&SUBMIT=Submit+request&option=1&optq=0' ...
          '&lat=49&lon=0&unit=1&plotom=0'];

% SETTINGS
stepsize = step;

% CREATING DATE STRINGS
dates = cell(2, 1);
for i = 1:2
    dates{i} = sprintf('&an%i=%i&mois%i=%i&jour%i=%i&h%i=%i&m%i=%i&s%i=%i', ...
        i, year(tspan(i)), ...
        i, month(tspan(i)), ...
        i, day(tspan(i)), ...
        i, hour(tspan(i)), ...
        i, minute(tspan(i)), ...
        i, round(second(tspan(i))));
end

% URL CREATING
url = [url, sprintf('step=%i', stepsize), dates{:}, params];

% READING URL
tic;
fprintf('Querying EOC...\n')
page = strsplit(webread(url), '\n').';
dt = toc;
fprintf('Query Time: %0.4f seconds\n', dt);

% EXTRACTING DATA
headers = convertCharsToStrings(strtrim(strsplit(page{2}, '/')));
dataStr = convertCharsToStrings(page(3:end-1));
data = zeros(length(dataStr), length(headers));
for i = 1:length(dataStr)
    data(i, :) = str2num(dataStr(i));
end

% FORMATTING DATA
s = struct;
for i = 1:length(dataStr)
    s(i).MJD = data(i, 1);
    s(i).Date = datetime(data(i, 2:7));
    s(i).R = reshape(data(i, 8:16), [3 3]);
    s(i).Quat = data(i, 17:end).';
end

end