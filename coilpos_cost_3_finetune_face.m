function cost = coilpos_cost_3_finetune_face(x, matsimnibs_orig, U, thetax_1, thetay_1, coil_center_2, target_skin, coilpos_base, coilpos_face, skin_msh, tmp_target_norm, coil_to_scalp_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Fine-tune coil position and orientation to minimize coil-target distance                  %
% Modified to use OPCODE for collision detection                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetax_3 = x(1); % Clockwise = positive
thetay_3 = x(2); % Clockwise = positive
zshift_3 = x(3);

rotMatx_3 = [1 0 0;0 cos(thetax_1+thetax_3),-sin(thetax_1+thetax_3);0,sin(thetax_1+thetax_3),cos(thetax_1+thetax_3)];
rotMaty_3 = [cos(thetay_1+thetay_3),0,sin(thetay_1+thetay_3);0,1,0;-sin(thetay_1+thetay_3),0,cos(thetay_1+thetay_3)];

matsimnibs_3 = matsimnibs_orig;
matsimnibs_3(1:3,1:3) = U*rotMatx_3*rotMaty_3*matsimnibs_orig(1:3,1:3);
coil_center_3 = coil_center_2+tmp_target_norm*zshift_3;
matsimnibs_3(1:3,4) = coil_center_3(:);

% Transform coil vertices using the transformation matrix
coilpos = (matsimnibs_3*[coilpos_base(:,1:3), ones(size(coilpos_base,1),1)]')';
coilpos = coilpos(:,1:3);

% Calculate minimum distance and check for collision using OPCODE
[min_distance, collision_flag] = get_coil_skin_distance(coilpos, coilpos_face, skin_msh, coil_to_scalp_distance);

% Calculate coil direction
coil_dir = (matsimnibs_3(1:3,1:3)*[0 0 1]')';
target_dir = target_skin-coil_center_3;
coil_angle = max(min(dot(coil_dir,target_dir)/(norm(coil_dir)*norm(target_dir)),1),-1);

% Calculate final cost
cost = exp(zshift_3+coil_to_scalp_distance) + exp((1-coil_angle)*10) ... 
    + (abs(min_distance-coil_to_scalp_distance)+1)^4 +(collision_flag+1)^4;
end