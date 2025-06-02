function angle = angle_between(vector1, vector2)
% ANGLE_BETWEEN Calculates the angle between two 3D vectors
%
% Syntax:
%   angle = angle_between(vector1, vector2)
%
% Inputs:
%   vector1 - First 3D vector (can be row or column)
%   vector2 - Second 3D vector (can be row or column)
%
% Outputs:
%   angle - Angle between the vectors in degrees
%
% Description:
%   This function calculates the angle between two 3D vectors using the
%   dot product formula. The result is returned in degrees.

% Ensure vectors are column vectors
vector1 = vector1(:);
vector2 = vector2(:);

% Normalize the vectors
unit_vector1 = vector1 / norm(vector1);
unit_vector2 = vector2 / norm(vector2);

% Calculate the dot product
dot_product = dot(unit_vector1, unit_vector2);

% Clamp the dot product to [-1, 1] to avoid numerical errors
dot_product = max(min(dot_product, 1), -1);

% Calculate the angle in radians
angle_rad = acos(dot_product);

% Convert to degrees
angle = rad2deg(angle_rad);

end