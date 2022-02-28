function bd = BodyData(id, opts)

% DOCS: https://ssd-api.jpl.nasa.gov/doc/horizons.html

arguments
    id 
    opts.Debug {mustBeNumericOrLogical} = false
    opts.CenterBody {mustBeInteger} = 0
end
%% INPUT HANDLING
if ismember(id, 1:9)
    id = 100*id + 99;
elseif ~isnumeric(id)
    if isstring(id); id = char(id); end
    idx = strcmpi({'sun', 'mercury', 'venus', 'earth', 'mars', ...
        'jupiter', 'saturn', 'uranus', 'neptune', 'pluto'}, id);
    idList = [10, 199:100:999];
    id = idList(idx);
    if isempty(idx)
        error(['Input string does not match list of bodies, ', ...
            'use body integer id''s.'])
    end
end

%% QUERYING HORIZONS
% INPUTS
dates(1) = datetime('today', 'Format', "uuuu-MMM-dd HH:mm:ss");
dates(2) = datetime(datenum(dates(1))-5*365, 'ConvertFrom', 'datenum');
dates = cellstr(datestr(dates, 'yyyy-mmm-dd HH:MM:SS'));

% CONSTRUCTING URL
url = sprintf(...
    ['https://ssd.jpl.nasa.gov/api/horizons.api?format=json', ...
    '&COMMAND=''%i''&OBJ_DATA=''YES''', ...
    '&EPHEM_TYPE=''ELEMENTS''&CENTER=''@%i''', ...
    '&START_TIME=''%s''&STOP_TIME=''%s''&STEP_SIZE=''1mo''', ...
    '&CSV_FORMAT=''YES'''], ...
    id, opts.CenterBody, dates{2}, dates{1});

% READING WEBSITE
tic
horizons = webread(url);
time = toc;
if opts.Debug
    horizons.result
end
txt = strsplit(horizons.result, '\n');

%% INITIALIZING DATA OUTPUT
bd = struct;
bd.Name = findBodyName(txt);

%% HANDLING CASES
switch id
    case 399  % Earth
        % DATA
        bd.Day = findNumeric(txt, 'Mean sidereal day', 1)*60;
        bd.Day_Solar = findNumeric(txt, 'Mean solar day 2000', 1);
        bd.J2 = findNumeric(txt, 'J2', 1);
        bd.MeanSMA = mean(findDataTable(txt, "A"));
        bd.MeanECC = mean(findDataTable(txt, "EC"));
        bd.MeanINC = mean(findDataTable(txt, "IN"));
        bd.Mu = findNumeric(txt, 'GM', 1);
        bd.Period = findNumeric(txt, 'Sidereal orb period', 2)*365;
        bd.Radius = findNumeric(txt, 'Equ. radius', 1);
        bd.SolarConstant = findNumeric(txt, 'Solar Constant', 1);
        
        % UNITS
        bd.Units.Day = 'sec';
        bd.Units.Day_Solar = 'sec';
        bd.Units.MeanSMA = 'km';
        bd.Units.MeanINC = 'deg';
        bd.Units.Mu = 'km3/s2';
        bd.Units.Period = 'days';
        bd.Units.Radius = 'km';
        bd.Units.SolarConstant = 'W/m2';
        
    otherwise
        % DATA
        searchStr = {'Sidereal rot. period', 'Mean solar day'};
        bd.Day = findNumeric(txt, searchStr, 1)*60;
        bd.J2 = findNumeric(txt, 'J2', 1);
        bd.MeanSMA = mean(findDataTable(txt, "A"));
        bd.MeanECC = mean(findDataTable(txt, "EC"));
        bd.MeanINC = mean(findDataTable(txt, "IN"));
        bd.Mu = findNumeric(txt, 'GM', 1);
        searchStr = {'Sidereal orb. per', 'Sidereal orb per', 'orbit period'};
        bd.Period = findNumeric(txt, searchStr, 1)*365;
        if isnan(bd.Period)
            bd.Period = mean(findDataTable(txt, "PR"))/86400;
        end
        searchStr = {'Vol. Mean Radius', 'Mean radius'};
        bd.Radius = findNumeric(txt, searchStr, 1);
        bd.SolarConstant = findNumeric(txt, 'Solar Constant', 1);
        
        % UNITS
        bd.Units.Day = 'sec';
        bd.Units.Day_Solar = 'sec';
        bd.Units.MeanSMA = 'km';
        bd.Units.MeanINC = 'deg';
        bd.Units.Mu = 'km3/s2';
        bd.Units.Period = 'days';
        bd.Units.Radius = 'km';
        bd.Units.SolarConstant = 'W/m2';
        
        
end
bd.QueryTime = time;

end

%% HELPER FUNCTIONS
function num = findNumeric(data, searchStr, ~)
    % CUTTING FROM DATA
    if iscell(searchStr)
        for i = 1:length(searchStr)
            idx = find(contains(data, searchStr{i}, "IgnoreCase", true), 1, 'first');
            if ~isempty(idx); searchStr = searchStr{i}; break; end
        end
    else
        idx = find(contains(data, searchStr, "IgnoreCase", true), 1, 'first');
    end
    if isempty(idx); num = NaN; return; end
    str = data{idx};
    
    % FINDING COLUMN
    dataIdx = strfind(lower(str), lower(searchStr));
    if length(dataIdx) > 1; dataIdx = dataIdx(1); end
    eqIdx = strfind(str, '=');
    if all(dataIdx < eqIdx)
        pos = 1;
    else
        pos = 2;
    end
    
    % FINDING DATA POSITION
    [s, e] = regexpi(str, '=\s+\d+([.]\d+(?=([\s\D]|$))|(?=\D))');
    if isempty(s); num = NaN; return; end
    num = str(s(pos):e(pos));
    
    % PULLING VALUE
    s = regexpi(num, '\d', "once");
    num = str2double(num(s:end));
end

function name = findBodyName(txt)
    % FINDING ROW
    idx = find(contains(txt, 'Revised: ', 'IgnoreCase', true), 1, 'first');
    if isempty(idx); name = 'could not find'; return; end
    name = txt{idx};
    
    % PULLING NAME
    [s, e] = regexpi(name, '\d+, \d+\s+\w*\s');
    if isempty(s); name = 'could not find'; return; end
    name = name(s:e);
    
    % CLEANING
    s = regexpi(name, '\s[a-z]+\0s');
    if isempty(s); name = 'could not find'; return; end
    name = name(s+1:end-1);
end

function val = findDataTable(txt, colname)
    % FINDING ROWS
    s = find(contains(txt, '$$SOE', 'IgnoreCase', true), 1, 'first')+1;
    e = find(contains(txt, '$$EOE', 'IgnoreCase', true), 1, 'first')-1;
    headers = textscan(txt{s-3}, '%s', 'Delimiter', ',');
    idx = strcmpi(headers{1}, colname);
    data = txt(s:e);
    val = zeros(length(data), length(colname));
    for i = 1:length(data)
        line = textscan(data{i}, '%s', 'Delimiter', ',');
        val(i) = str2double(line{1}{idx});
    end
end