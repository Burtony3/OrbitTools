function angle = u2vAng(u, v, units, type)

% ARGUMENT VALIDATION
arguments
    u {mustBeNumeric}
    v {mustBeNumeric}
    units {mustBeTextScalar} = 'degrees'
    type {mustBeTextScalar, mustBeMember(type, ...
        ["clockwise", "inside", "counterclockwise"])} = 'clockwise'
end

% Units
if strncmpi(units, 'degrees', 1)
    units = 180/pi;
else
    units = 1;
end

% Angle Type
if strcmpi(type, 'inside')
    func = @(angle) -angle;
else
    func = @(angle) 2*pi - angle;
end
    

% Handling Column Vecs
if all(size(u) == [3 1])
    u = u(:).';
    v = v(:).';
end

% CALCULATING ANLGES
angle = zeros(size(u, 1), 1);
for i = 1:size(u, 1)
    % ENFORCING COLUMN VECS
    u_ = u(i, :).';
    v_ = v(i, :).';

    % FINDING VECTOR PRODUCTS
    c = cross(u_, v_);
    d = dot(u_, v_);

    % FINDING ANGLE
    angle(i) = atan2(norm(c), d);

    % ACCOUNTING FOR QUADRENT
    if c(3) < 0
        angle(i) = func(angle(i));
    end
end

% OUTPUTTING
angle = angle*units;

end