function R = rot(type, angles, units)


%% INPUT VALIDATION
% ERROR CATCHING
type = char(type);
if length(char(type)) ~= length(angles)
    error('Rotation type length does not match number of angles')
end
if ~all(type == '1' | type == '2' | type == '3')
    error('Valid rotation types are must be in form of ''321'', or some permutation')
end

if strncmpi(units, 'd', 1)
    units = pi/180;
else
    units = 1;
end

%% CREATING MATRIX
% HELPFUL VALUES
nangles = length(angles);

% CREATING MATRIX
R = eye(3);
for i = 1:nangles
    R = R*rot_sub(angles(i)*units, type(i));
end

%% Rotation Function
function R = rot_sub(angle, num)
    switch num
        case '1'
            R = [ 1  0            0          ;
                  0  cos(angle)  sin(angle);
                  0 -sin(angle)  cos(angle) ];
        case '2'
            R = [ cos(angle)  0 -sin(angle); 
                  0           1  0          ;
                  sin(angle)  0  cos(angle) ];
            
        case '3'
            R = [ cos(angle)  sin(angle)  0;
                 -sin(angle)  cos(angle)  0;
                  0            0            1 ];
            
    end
end

end