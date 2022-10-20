function varargout = coe2rv(coe, mu)
    a = coe(1);
    e = coe(2);
    i = coe(3)*pi/180;
    Omega = coe(4)*pi/180;
    omega = coe(5)*pi/180;
    TA = coe(6)*pi/180;
    
    r = a*(1-e^2) / (1 + e*cos(TA));
    p = a*(1-e^2);
    h = sqrt(mu * p);
    
    rvec = r * [...
        cos(Omega)*cos(omega+TA) - sin(Omega)*sin(omega+TA)*cos(i); 
        sin(Omega)*cos(omega+TA) + cos(Omega)*sin(omega+TA)*cos(i);
        sin(i)*sin(omega+TA)];
    
    vvec = [...
        ((rvec(1)*h*e)/(r*p))*sin(TA) - (h/r)*(cos(Omega)*sin(omega+TA) + sin(Omega)*cos(omega+TA)*cos(i));
        ((rvec(2)*h*e)/(r*p))*sin(TA) - (h/r)*(sin(Omega)*sin(omega+TA) - cos(Omega)*cos(omega+TA)*cos(i));
        ((rvec(3)*h*e)/(r*p))*sin(TA) + (h/r)*(sin(i)*cos(omega+TA))];

    if nargout == 1
        varargout{1} = [rvec; vvec];
    else
        varargout{1} = rvec;
        varargout{2} = vvec;
    end

end