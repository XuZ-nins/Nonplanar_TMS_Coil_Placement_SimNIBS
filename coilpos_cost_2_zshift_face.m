function cost = coilpos_cost_2_zshift_face(x, matsimnibs_1, coilpos_base, coilpos_face, skin_msh, target_skin, tmp_target_norm, coil_to_scalp_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Shift coil position such that the coil-scalp distance is equal to the preset value        %
% Modified to use OPCODE for collision detection                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zshift_2 = x(1);

matsimnibs_2 = matsimnibs_1;
coil_center_2 = target_skin + tmp_target_norm * (zshift_2 + coil_to_scalp_distance);
matsimnibs_2(1:3,4) = coil_center_2(:); % Ensure no intrusion at first

% Transform coil vertices using the transformation matrix
coilpos = (matsimnibs_2*[coilpos_base(:,1:3), ones(size(coilpos_base,1),1)]')';
coilpos = coilpos(:,1:3);

% Calculate minimum distance and check for collision using OPCODE
[min_distance, collision_flag] = get_coil_skin_distance(coilpos, coilpos_face, skin_msh, coil_to_scalp_distance);

% Calculate final cost
cost = (abs(min_distance-coil_to_scalp_distance)+1)^4 + (collision_flag+1)^4;
end