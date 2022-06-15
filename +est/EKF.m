function [xkhat, Pk] = EKF(f, x, P, zk, h, Q, R)

% GENERATING STATE TRANSITION MATRIX
Fk = opt.finiteDiff(f, {x}, 1, h = 1e-6, szOut = 6).';

% GENERATING PREDICTED STATE AND COVARIANCE
xphat = f(x);
Pp = Fk*P*Fk.' + Q;

% GENERATING OBSERVATION MATRIX
Hk = opt.finiteDiff(h, {xphat}, 1, h = 1e-6, szOut = 2).';

% KALMAN GAIN
yk = zk - h(xphat);             % Innovation
Sk = Hk*Pp*Hk.' + R;            % Innovation Covariance
Kk = Pp*Hk.'*inv(Sk);           % Kalman Gain

% STATE UPDATE EQUATIONS
xkhat = xphat + Kk*yk;
Pk = (eye(length(x)) - Kk*Hk)*Pp;

end