function Yf = KepUV( x0,deltaUV)
% Ryan P. Russell  3-23-2016
% INPUTS:
% x0        initial position and velocity
% deltaUV   change in the universal varialbe (independent variable of the propagation)
% [REMOVED] orderCase what order derivative 0,1 or 2

% OUTPUTS:
% Yf        final position, velocity, and DeltaTime
% [REMOVED] dYf       (only populates if orderCase>0) partial of Y wrt initial position, initial velocity, and deltaUV
% [REMOVED] d2Yf      (only populates if orderCase>1) second partial of Y wrt initial position, initial velocity, and deltaUV
% [REMOVED] capZ      value of z from BMW at output (negative hyperbola, positive ellipse, near zero is near parabola
% [REMOVED] errFlag   0 if all ok, -100 if invalid OrderCase, -1 if too hyperbolic such that sinh and cosh will overflow: if arg bigger than ~sqrt(504100)

%you need to initialize theses
Yf=zeros(7,1);
p=zeros(7,1);
%=================================================
%below are params, only need to be set once
epsswitch=0.01d0;  %the maple code is designed to work with 0.01; meaning the terms in the C and S series are hard coded to be long enough so that when z=0.01, the result is the same to 16 digits when using the full expressions.
NEGepsswitch=-epsswitch;
%onebyQfact(1:9)=[1.d0/1.d0, 1.d0/2.d0,1.d0/6.d0,1.d0/24.d0,1.d0/120.d0,1.d0/720.d0,1.d0/5040.d0,1.d0/40320.d0,1.d0/362880.d0];
onebyQfact(1:9)=[1,0.5,0.166666666666667,0.0416666666666667,0.00833333333333333,0.00138888888888889,0.000198412698412698,2.48015873015873e-05,2.75573192239859e-06];
%=================================================
x=deltaUV;


xsq=x * x;
v0sq=x0(4)*x0(4)+x0(5)*x0(5)+x0(6)*x0(6);
r0sq=x0(1)*x0(1)+x0(2)*x0(2)+x0(3)*x0(3);
r0 = sqrt(r0sq);
onebyr0=1.d0/r0;
capZ = xsq * (-v0sq * r0 + 2.d0)  * onebyr0;

%check for overflow condition (sinh and cosh will overflow if arg bigger than ~sqrt(504100) too big
if abs(capZ)>50400.d0
    Yf=0.d0;
    return
end


if capZ>epsswitch
    zsign=1;
elseif (capZ<NEGepsswitch)
    zsign=-1;
else
    zsign=0;
end

if zsign ==   -1
    t1 = sqrt(-capZ);
    t2 = cosh(t1);
    t3 = 1 - t2;
    t5 = -r0 * v0sq + 2;
    t6 = 1 / t5;
    t7 = t6 * t3;
    t8 = 1 - t7;
    t13 = x0(1) * x0(4) + x0(2) * x0(5) + x0(3) * x0(6);
    t15 = r0 * t6;
    t16 = t15 * t3 * t13;
    t18 = sinh(t1);
    t22 = -0.1e1 / t1 / capZ * (t18 - t1);
    t26 = -onebyr0 * t5 * xsq * t22 + 1;
    t27 = t26 * x * r0;
    t28 = t16 + t27;
    t42 = 0.1e1 / (t26 * x * t13 + t2 * r0 + r0 * t7);
    t43 = t42 * onebyr0;
    t44 = -t26 * x;
    t49 = -t15 * t3 * t42 + 1;
    p(1) = x0(4) * t28 + x0(1) * t8;
    p(2) = x0(5) * t28 + x0(2) * t8;
    p(3) = x0(6) * t28 + x0(3) * t8;
    p(4) = x0(1) * t44 * t43 + x0(4) * t49;
    p(5) = x0(2) * t44 * t43 + x0(5) * t49;
    p(6) = x0(3) * t44 * t43 + x0(6) * t49;
    p(7) = t22 * xsq * x + t16 + t27;
    Yf(1:7)=p(1:7);
          
elseif zsign == 0
    t3 = -r0 * v0sq + 2;
    t4 = t3 * xsq;
    t7 = xsq ^ 2;
    t8 = t3 ^ 2;
    t9 = t8 * t7;
    t10 = 0.1e1 / r0sq;
    t13 = x ^ 2;
    t14 = t13 ^ 2;
    t17 = t8 * t3 * t14 * t13;
    t18 = t10 * onebyr0;
    t21 = 0.1e1 / 0.2e1 - onebyr0 * onebyQfact(4) * t4 + onebyQfact(6) * t10 * t9 - onebyQfact(8) * t18 * t17;
    t23 = -t21 * onebyr0 * xsq + 1;
    t28 = x0(1) * x0(4) + x0(2) * x0(5) + x0(3) * x0(6);
    t30 = t21 * xsq * t28;
    t38 = 0.1e1 / 0.6e1 - onebyr0 * onebyQfact(5) * t4 + onebyQfact(7) * t10 * t9 - onebyQfact(9) * t18 * t17;
    t40 = onebyr0 * t3;
    t42 = -t40 * xsq * t38 + 1;
    t43 = t42 * x * r0;
    t44 = t30 + t43;
    t53 = t21 * xsq;
    t60 = 0.1e1 / (t53 + t42 * x * t28 + (-t40 * t53 + 1) * r0);
    t61 = t60 * onebyr0;
    t62 = -t42 * x;
    t67 = -t21 * t60 * xsq + 1;
    p(1) = x0(1) * t23 + x0(4) * t44;
    p(2) = x0(2) * t23 + x0(5) * t44;
    p(3) = x0(3) * t23 + x0(6) * t44;
    p(4) = x0(1) * t62 * t61 + x0(4) * t67;
    p(5) = x0(2) * t62 * t61 + x0(5) * t67;
    p(6) = x0(3) * t62 * t61 + x0(6) * t67;
    p(7) = t38 * x * xsq + t30 + t43;
    Yf(1:7)=p(1:7);

elseif zsign == 1
    t1 = sqrt(capZ);
    t2 = cos(t1);
    t3 = 1 - t2;
    t5 = -r0 * v0sq + 2;
    t6 = 1 / t5;
    t7 = t6 * t3;
    t8 = 1 - t7;
    t13 = x0(1) * x0(4) + x0(2) * x0(5) + x0(3) * x0(6);
    t15 = r0 * t6;
    t16 = t15 * t3 * t13;
    t18 = sin(t1);
    t22 = 0.1e1 / t1 / capZ * (t1 - t18);
    t26 = -onebyr0 * t5 * xsq * t22 + 1;
    t27 = t26 * x * r0;
    t28 = t16 + t27;
    t42 = 0.1e1 / (t26 * x * t13 + t2 * r0 + r0 * t7);
    t43 = t42 * onebyr0;
    t44 = -t26 * x;
    t49 = -t15 * t3 * t42 + 1;
    p(1) = x0(4) * t28 + x0(1) * t8;
    p(2) = x0(5) * t28 + x0(2) * t8;
    p(3) = x0(6) * t28 + x0(3) * t8;
    p(4) = x0(1) * t44 * t43 + x0(4) * t49;
    p(5) = x0(2) * t44 * t43 + x0(5) * t49;
    p(6) = x0(3) * t44 * t43 + x0(6) * t49;
    p(7) = t22 * xsq * x + t16 + t27;
    Yf(1:7)=p(1:7);

end
% Yf = norm(Yf);
% Yf = sqrt(sum(Yf.^2));
end


