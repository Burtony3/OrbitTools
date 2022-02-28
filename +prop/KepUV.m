function [Yf,dYf,d2Yf,capZ,errFlag] = propKeplerUV( x0,deltaUV,orderCase,loadLib,closeLib)
% this calls a fortran dynamic library 
% both 'kepUV.dll' and 'kepUV.h' must be in your matlab path
% Ryan P. Russell  3-23-2016 
% INPUTS:
% x0        initial position and velocity
% deltaUV   change in the universal varialbe (independent variable of the propagation)
% orderCase what order derivative 0,1 or 2
% loadLib   true or false (open library first)   DOES NOT COMPUTE TRAJ, JUST OPENS LIB
% closeLib  true of false (close library when done)  DOES NOT COMPUTE TRAJ, JUST CLOSES LIB

% OUTPUTS:
% Yf        final position, velocity, and DeltaTime
% dYf       (only populates if orderCase>0) partial of Y wrt initial position, initial velocity, and deltaUV
% d2Yf      (only populates if orderCase>1) second partial of Y wrt initial position, initial velocity, and deltaUV
% capZ      value of z from BMW at output (negative hyperbola, positive ellipse, near zero is near parabola
% errFlag   0 if all ok, -100 if invalid OrderCase, -1 if too hyperbolic such that sinh and cosh will overflow: if arg bigger than ~sqrt(504100)

%make these global so that you only need to initialize them once, when
% global Yf dYf d2Yf capZ errFlag

% Step 1: Load the shared library
% The intrinsic MATLAB function 'loadlibrary' takes in 2 arguments:
% First, the shared library *.dll file (including directory path, if
% needed). Second, the shared library *.h header file, in which prototypes
% of the routines in the *.dll are given, C-style (including directory
% path, if necessary).
% The library '*.dll' is loaded as simply '*'.



%q=length(Yf);
    
Yf=NaN(7,1);
dYf=NaN(7,7);
d2Yf=NaN(7,7,7);
capZ=NaN;
errFlag=NaN;

if(loadLib)
    
    Yf=NaN(7,1);
    dYf=NaN(7,7);
    d2Yf=NaN(7,7,7);
    capZ=NaN;
    errFlag=NaN;
    
    
    %*********EDIT ARGUMENTS TO REFLECT DIRECTORY STRUCTURE
    if ~libisloaded('kepUV') % if library is not loaded, load it
        loadlibrary('kepUV.dll','kepUV.h'); % load on client
    else % unload it, then load it
        unloadlibrary('kepUV') % unload on client
        loadlibrary('kepUV.dll','kepUV.h'); % load on client
    end
    
%     disp('kepUV.dll library loaded')
    
    % The instrinsic MATLAB function 'libfunctions' can be used to view the
    % functions within the library that are now user-callable.
    
    %uncomment if you want to list
    %list = libfunctions('kepUV', '-full');
    
    
end

if(not(loadLib)&not(closeLib))
    if(orderCase>2|orderCase<0)
        disp('error, invalid orderCase in propKeplerUV')
        errFlag=-100
        return
    end
    % call routine
    % first string arg of calllib is name of library, second is name of routine
    [~,~,Yf,dYf,d2Yfo,capZ,~,errFlag] = calllib('kepUV','PartialsKepDeltaE', x0,deltaUV,Yf,dYf,d2Yf,capZ,orderCase,errFlag); % inputs
    d2Yf=reshape(d2Yfo,[7,7,7]);  %some reason matlab thinks its a 7x49 matrix not a tensor?
    if(errFlag==-1)
        disp('warning: propKeplerUV cant propagate in double precision because solution is too hyperbolic')
    end
    
end

if(closeLib)
    %The intrinsic MATLAB function 'unloadlibrary' unloads the shared library.
    unloadlibrary('kepUV') % unload on client
%     disp('kepUV.dll library unloaded')
end

end

