function cost = coilpos_cost_1_xyangle(x,matsimnibs_orig,coil_center_1,U,shp,coilpos_base,skin_norm,coil_to_scalp_distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Find initial coil orientation that likely minimizes coil-target distance                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetax_1 = x(1); % Clockwise = positive
thetay_1 = x(2); % Clockwise = positive

rotMatx_1 = [1 0 0;0 cos(thetax_1),-sin(thetax_1);0,sin(thetax_1),cos(thetax_1)];
rotMaty_1 = [cos(thetay_1),0,sin(thetay_1);0,1,0;-sin(thetay_1),0,cos(thetay_1)];

matsimnibs_1(1:3,4) = coil_center_1(:);
matsimnibs_1(1:3,1:3) = U*rotMatx_1*rotMaty_1*matsimnibs_orig(1:3,1:3);

coilpos = (matsimnibs_1*coilpos_base')';
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

cost = -minD + (inShapeCost+1)^4 - 1;

end