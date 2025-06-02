function cost = coilpos_cost_3_finetune(x,matsimnibs_orig,U,thetax_1,thetay_1,coil_center_2,target_skin,shp,coilpos_base,tmp_target_norm,skin_norm,coil_to_scalp_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Fine-tune coil position and orientation to minimize coil-target distance                  %
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
coilpos = (matsimnibs_3*coilpos_base')';
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

coil_dir = (matsimnibs_3(1:3,1:3)*[0 0 1]')';
target_dir = target_skin-coil_center_3;
coil_angle = max(min(dot(coil_dir,target_dir)/(norm(coil_dir)*norm(target_dir)),1),-1);

% cost = exp(norm(target_ROI-coil_center_3-zrange/2)/100) + (coil_angle>=0)*(1-coil_angle)/2 + (coil_angle<0)*(1-coil_angle)*100 ...
%     + (minD>coil_to_scalp_distance)*(minD-coil_to_scalp_distance+1)^2 ...
%     + (minD<=coil_to_scalp_distance)*(coil_to_scalp_distance-minD+1)^2 -1 ...
%     + (inShapeCost+1)^4*1 - 1;

% +(1-coil_angle)/2

cost = exp(zshift_3+coil_to_scalp_distance) + exp((1-coil_angle)*10) ... % (coil_angle>=0)*(1-coil_angle)/2 + (coil_angle<0)*(1-coil_angle)*100 ...
    + (abs(minD-coil_to_scalp_distance)+1)^4 - 1 + (inShapeCost+1)^4 - 1;

end
