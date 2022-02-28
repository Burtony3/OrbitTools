function coe = rv2coe(rvec, vvec, opts)
    arguments
        rvec {mustBeNumeric} 
        vvec {mustBeNumeric} 
        opts.Mu {mustBeNumeric, mustBePositive} = 1
        opts.Anomaly {mustBeText, mustBeMember(opts.Anomaly, {'true', 'eccentric', 'mean'})} = 'true'
        opts.Units {mustBeText, mustBeMember(opts.Units, {'deg', 'rad'})} = 'deg'
        opts.Dimension {mustBeInteger, mustBeInRange(opts.Dimension, 2, 3)} = 3
        opts.PrintElements {mustBeNumericOrLogical} = false
    end
    
    % HANDLING OPTIONALS
    anomaly = find(strcmpi(opts.Anomaly, {'true', 'eccentric', 'mean'}));  % Anomaly Type
    if strcmpi(opts.Units, 'deg'); units = 180/pi; else; units = 1; end    % Output Units
    
    % CORRECTING FOR INPUT SIZE
    sz = size(rvec);
    if ~any(sz == opts.Dimension)
        error('Input size does not match dimension\n Size = [%i %i]\n Expected Dimension = %i',...
            sz(1), sz(2), opts.Dimension)
    end
    if ~all(sz == size(vvec))
        error('Vector input dimension mismatch between velocity and radius'); 
    end
    if sz(1) ~= opts.Dimension
        rvec = rvec.';
        vvec = vvec.';
    end
    
    % PREFILLING
    height = sz(sz ~= opts.Dimension);
    coe = zeros(height, opts.Dimension*2);
    
    % LOOPING
    for j = 1:height 
        coe(j, :) = subfunc(rvec(:, j), vvec(:, j), opts.Mu, anomaly, units, opts.Dimension);
    end
    
    if opts.PrintElements
        anomSym = [957, 69, 77];
        fprintf('<strong>coe = [a, e, i, %c, %c, %c]</strong>\n', 937, 969, anomSym(anomaly))
    end
    
    function coe = subfunc(rvec, vvec, mu, anomaly, units, dim)
        % FINDING NORMS
        r = norm(rvec);
        v = norm(vvec);
        
        % HANDLING IN PLANE SINGULARITIES
        if rvec(3) == 0 && vvec(3) == 0
            rvec(3) = 1e-12;
        end
        
        % HANDLING DIMENSIONS
        if dim == 2
            rvec = [rvec; 1e-12];
            vvec = [vvec; 0];
        end

%         % SEMI-MAJOR AXIS
%         a = 1/(2/r - (v^2)/mu);

        % SPECIFIC ANGULAR MOMENTUM
        hvec = cross(rvec, vvec);
        h = norm(hvec);

        % ECCENTRICITY VECTOR
        evec = cross(vvec, hvec)/mu - rvec/r;
        e = norm(evec);
        
        % SEMI-MAJOR AXIS Ver. 2
        epsilon = 0.5*v^2 - mu/r;
        a = mu / (2*epsilon);
        
        % INCLINATION
        i = acos(hvec(3)/h);

        % LINE OF NODES
        nvec = [-hvec(2); hvec(1); 0];
        n = norm(nvec);

        % LONGITUDE OF ASCENDING NODE
        if nvec(2) >= 0
            Omega = acos(nvec(1)/n);
        else
            Omega = 2*pi - acos(nvec(1)/n);
        end

        % ARGUMENT OF PERIAPSIS
        if evec(3) >= 0
            omega = acos(dot(nvec, evec)/(n*e));
        else
            omega = 2*pi - acos(dot(nvec, evec)/(n*e));
        end

        % TRUE ANOMALY
        tmp = min(1, dot(evec, rvec)/(e*r));
        if dot(rvec, vvec) > 0
            theta = acos(tmp);
        else
            theta = 2*pi - acos(tmp);
        end
        angle_out = theta;

        % ECCENTRIC ANOMALY
        if anomaly > 1
            if e < 2
                E = 2*atan(tan(theta/2)/sqrt(1 + e/(1 - e)));
            else
                E = 2*atanh(tan(theta/2)/sqrt(1 - e/(1 - e)));
            end
            angle_out = E;
        end

        % MEAN ANOMALY
        if anomaly > 2
            if e < 1
                M = E - e*sin(E);
            else
                M = e*sinh(E) - E;
            end
            angle_out = M;
        end

        % OUTPUTTING
        if dim == 3
            coe = [a, e, i*units, Omega*units, omega*units, angle_out*units].';
        elseif dim == 2
            coe = [a, e, omega*units, angle_out*units].';
        end
    end
end