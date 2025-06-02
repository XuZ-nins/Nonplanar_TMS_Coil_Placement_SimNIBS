function cost = coilpos_cost_1_xyangle_face(x, matsimnibs_orig, coil_center_1, U, coilpos_base, coilpos_face, skin_msh, target_ROI, coil_to_scalp_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Find initial coil orientation that likely minimizes coil-target distance                  %
% Modified to use OPCODE for collision detection                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetax_1 = x(1); % Clockwise = positive
thetay_1 = x(2); % Clockwise = positive

rotMatx_1 = [1 0 0;0 cos(thetax_1),-sin(thetax_1);0,sin(thetax_1),cos(thetax_1)];
rotMaty_1 = [cos(thetay_1),0,sin(thetay_1);0,1,0;-sin(thetay_1),0,cos(thetay_1)];

matsimnibs_1 = eye(4);
matsimnibs_1(1:3,4) = coil_center_1(:);
matsimnibs_1(1:3,1:3) = U*rotMatx_1*rotMaty_1*matsimnibs_orig(1:3,1:3);

% Transform coil vertices using the transformation matrix
coilpos = (matsimnibs_1*[coilpos_base(:,1:3), ones(size(coilpos_base,1),1)]')';
coilpos = coilpos(:,1:3);

coil_vec_1 = transform_normal(matsimnibs_1, [0,0,1]);
ROI_vec_1 = target_ROI-coil_center_1; ROI_vec_1 = ROI_vec_1/norm(ROI_vec_1);

% Check for infeasible solution 
ROI_coil_angle = angle_between(ROI_vec_1,coil_vec_1); % Coil focus deviates too much from ROI

% Calculate minimum distance and check for collision using OPCODE
[min_distance, collision_flag] = get_coil_skin_distance(coilpos, coilpos_face, skin_msh, coil_to_scalp_distance);

% Calculate final cost
cost = -min_distance + ROI_coil_angle_penalty(ROI_coil_angle)*10 + (collision_flag+1)^4;
end

function y = ROI_coil_angle_penalty(x)
% Set default parameters if not provided
a = 0.05;  % Base growth rate
b = 0.15;  % Accelerated growth factor
c = 15;    % Transition point

% Calculate sigmoid transition function (smooth step from 0 to 1)
sigmoid = 1 ./ (1 + exp(-0.3 * (x - c)));

% Calculate base exponential component
base_exp = exp(a * x);

% Calculate enhanced exponential growth that increases smoothly after x=c
enhanced_exp = exp(b * sigmoid .* (x - c) + a * x);

% Smoothly blend between base and enhanced exponential 
% to create continuous acceleration
y = base_exp .* (1 - sigmoid) + enhanced_exp .* sigmoid - 1;

end