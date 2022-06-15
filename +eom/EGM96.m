function dX = EGM96(~, X)
    
% CONVERTING STATE KM -> M & COL -> ROW
r = X(1:3).'*1000;

% OUTPUTTING EGM MODEL
[gx, gy, gz] = gravitysphericalharmonic(r, "EGM96", 20); % Input/output in m

% CREATING ACCELERATION VECTOR
dX = [X(4:6);
      [gx; gy; gz]/1000]; % Back to km/s^2

end