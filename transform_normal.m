function transformed_normal = transform_normal(transformation_matrix, normal_vector)
% TRANSFORM_NORMAL Transforms a normal vector using a 4x4 transformation matrix
%
% Syntax:
%   transformed_normal = transform_normal(transformation_matrix, normal_vector)
%
% Inputs:
%   transformation_matrix - 4x4 transformation matrix
%   normal_vector - 3x1 or 1x3 normal vector to be transformed
%
% Outputs:
%   transformed_normal - The transformed normal vector (normalized)
%
% Description:
%   This function transforms a normal vector using the correct method,
%   which is to multiply the normal by the transpose of the inverse of the
%   upper-left 3x3 submatrix of the transformation matrix.

% Ensure normal vector is a column vector
normal_vector = normal_vector(:);

% Extract the 3x3 rotation/scaling part from the transformation matrix
R = transformation_matrix(1:3, 1:3);

% Calculate the inverse of the 3x3 matrix
R_inv = inv(R);

% Calculate the transpose of the inverse
R_inv_transpose = R_inv';

% Apply the transformation to the normal vector
transformed_normal = R_inv_transpose * normal_vector;

% Normalize the result
transformed_normal = transformed_normal / norm(transformed_normal);

end
