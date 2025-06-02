function cost = coilpos_cost_4_closeup(x,matsimnibs_3,shp,coilpos_base,tmp_target_norm,skin_norm,coil_to_scalp_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Final optimization step to ensure the coil-scalp distance is equal to the preset value    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zshift_4 = x(1); % Clockwise = positive

matsimnibs_final = matsimnibs_3;
coil_center_4 = matsimnibs_3(1:3,4)'+tmp_target_norm*zshift_4; % Ensure no intrusion at first
matsimnibs_final(1:3,4) = coil_center_4(:);

coilpos = (matsimnibs_final*coilpos_base')';
[I,D] = nearestNeighbor(shp,coilpos(:,1:3));

% Calibrate distance due to loss of resolution
minDind = 1:length(D);
% minDind = D<min(D)+1;
minD_list = D(minDind);
minDcoil = coilpos(minDind,1:3);
minDPoint = shp.Points(I(minDind),:);

a = minDcoil-minDPoint;
b = skin_norm(I(minDind),:);
angle_list = sum(a.*b,2)./(vecnorm(a,2,2).*vecnorm(b,2,2));
valid_pair = (abs(angle_list)>sqrt(2)/2 | minD_list<coil_to_scalp_distance) & minD_list<coil_to_scalp_distance*3;
minD_true = minD_list.*angle_list;
minD_valid = minD_true(valid_pair);
minD = min(minD_valid(minD_valid>0));
if min(abs(minD_true(angle_list<0)))<coil_to_scalp_distance*2
    inShapeCost = sum(angle_list<0);
else
    inShapeCost = 0;
end

if isempty(minD)
    minD = min(minD_list);
end

cost = (abs(minD-coil_to_scalp_distance)+1)^4 - 1 + (inShapeCost+1)^4 - 1;

end