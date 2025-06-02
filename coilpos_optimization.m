function [matsimnibs_final,D_ROI,D_skin] = coilpos_optimization(shp,coilpos_base,U,thetaz,maxz0,zrange,...
    target_ROI,target_skin,skin_norm,tmp_target_norm,coil_to_scalp_distance,xangle_range,yangle_range,xyangle_tolerance,fine_tuning)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Main function for coil position optimization                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matsimnibs_orig = zeros(4,4);
matsimnibs_orig(4,4) = 1;
coil_center_1 = target_skin+tmp_target_norm*(max(zrange*2,30)+coil_to_scalp_distance); % Take coil afar at first
rotMatz_1 = [cos(thetaz),-sin(thetaz) 0;sin(thetaz),cos(thetaz),0;0,0,1];
matsimnibs_orig(1:3,1:3) = rotMatz_1;

% Optimization 1: Find best x/y-axis rotations such that the
% coil faces the scalp with roughly equal distances across both
% wings (by maximing the total coil-scalp distance)
f = @(x)coilpos_cost_1_xyangle(x,matsimnibs_orig,coil_center_1,U,shp,coilpos_base,skin_norm,coil_to_scalp_distance);
x0 = [0 0];
lb = [-xangle_range -yangle_range];
ub = [xangle_range yangle_range];
options = optimoptions(@fmincon,'MaxFunctionEvaluations',1e3,'StepTolerance',xyangle_tolerance, ...
    'Display','off'); % xyangle_tolerance tolerance
[x,fval1] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);

thetax_1 = x(1); % Clockwise = positive
thetay_1 = x(2); % Clockwise = positive
rotMatx_1 = [1 0 0;0 cos(thetax_1),-sin(thetax_1);0,sin(thetax_1),cos(thetax_1)];
rotMaty_1 = [cos(thetay_1),0,sin(thetay_1);0,1,0;-sin(thetay_1),0,cos(thetay_1)];
matsimnibs_1 = matsimnibs_orig;
matsimnibs_1(1:3,1:3) = U*rotMatx_1*rotMaty_1*matsimnibs_1(1:3,1:3);
matsimnibs_1(1:3,4) = coil_center_1';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Bypass Optimization 1
% matsimnibs_1(1:3,4) = coil_center_1';
% matsimnibs_1(1:3,1:3) = U*matsimnibs_orig(1:3,1:3);
% 
% coilpos = (matsimnibs_1*coilpos_base')';
% [I,D] = nearestNeighbor(shp,coilpos(:,1:3));
% 
% % Calibrate distance due to loss of resolution
% minDind = 1:length(D);
% % minDind = D<min(D)+1;
% minD_list = D(minDind);
% minDcoil = coilpos(minDind,1:3);
% minDPoint = shp.Points(I(minDind),:);
% 
% a = minDcoil-minDPoint;
% b = skin_norm(I(minDind),:);
% angle_list = sum(a.*b,2)./(vecnorm(a,2,2).*vecnorm(b,2,2));
% valid_pair = (abs(angle_list)>sqrt(2)/2 | minD_list<coil_to_scalp_distance) & minD_list<coil_to_scalp_distance*3;
% minD_true = minD_list.*angle_list;
% minD_valid = minD_true(valid_pair);
% D_skin = min(minD_valid(minD_valid>0));
% if min(abs(minD_true(angle_list<0)))<coil_to_scalp_distance*2
%     inShapeCost = sum(angle_list<0);
% else
%     inShapeCost = 0;
% end
% 
% if isempty(minD)
%     minD = min(minD_list);
% end
% 
% fval1 = -minD + (inShapeCost+1)^4 - 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false % fval1>0 % Infeasible
    coilpos = (matsimnibs_1*coilpos_base')';
    fprintf('Cost = Inf\n');
    matsimnibs_final = zeros(4,4);
    D_ROI = Inf;
    D_skin = Inf;
else    
    % Optimization 2: Find best shift along the normal vector of
    % the current skin vertex, so as to minimize the coil-scalp distance
    % 1-D optimization (fminbnd); no constraints necessary
    f = @(x)coilpos_cost_2_zshift(x,matsimnibs_1,shp,coilpos_base,target_skin,tmp_target_norm,skin_norm,coil_to_scalp_distance);
    options = optimset('MaxIter',100,'MaxFunEvals',100);
    lb = maxz0;
    ub = max(zrange*2,30)+maxz0; % Planar coils need to stay further...
    [x,fval2] = fminbnd(f,lb,ub,options);
    if fval2>10 && ~fine_tuning % Brute-force approach only when fine_tuning is disabled
        cost_all = zeros(100,1);
        zshift_all = linspace(lb,ub,100);
        for ii = 1:length(zshift_all)
            zshift_tmp = zshift_all(ii);
            cost_all(ii) = f(zshift_tmp);
        end
        [fval2,Iminc] = min(cost_all);
        x = zshift_all(Iminc);
    end
    zshift_2 = x(1);
    matsimnibs_2 = matsimnibs_orig;
    matsimnibs_2(1:3,1:3) = U*rotMatx_1*rotMaty_1*matsimnibs_orig(1:3,1:3);
    coil_center_2 = target_skin+tmp_target_norm*(zshift_2+coil_to_scalp_distance);
    matsimnibs_2(1:3,4) = coil_center_2(:); % Ensure no intrusion at first

    if fine_tuning
        % Optimization 3: Precise tuning of x/y rotations along with z shifts
        f = @(x)coilpos_cost_3_finetune(x,matsimnibs_orig,U,thetax_1,thetay_1,coil_center_2,target_skin,shp,coilpos_base,tmp_target_norm,skin_norm,coil_to_scalp_distance);
        x0 = [0 0 0];
        lb = [-xangle_range -yangle_range -coil_to_scalp_distance-zrange/4];
        ub = [xangle_range yangle_range coil_to_scalp_distance];
        options = optimoptions(@fmincon,'MaxFunctionEvaluations',1e3,'StepTolerance',1e-3, ...
            'Display','off'); % xyangle_tolerance tolerance
        [x,fval3] = fmincon(f,x0,[],[],[],[],lb,ub,[],options);
        thetax_3 = x(1); % Clockwise = positive
        thetay_3 = x(2); % Clockwise = positive
        zshift_3 = x(3);
        
        rotMatx_3 = [1 0 0;0 cos(thetax_1+thetax_3),-sin(thetax_1+thetax_3);0,sin(thetax_1+thetax_3),cos(thetax_1+thetax_3)];
        rotMaty_3 = [cos(thetay_1+thetay_3),0,sin(thetay_1+thetay_3);0,1,0;-sin(thetay_1+thetay_3),0,cos(thetay_1+thetay_3)];
        
        matsimnibs_3 = matsimnibs_orig;
        matsimnibs_3(1:3,1:3) = U*rotMatx_3*rotMaty_3*matsimnibs_orig(1:3,1:3);
        coil_center_3 = coil_center_2+tmp_target_norm*zshift_3;
        matsimnibs_3(1:3,4) = coil_center_3(:);

        % Optimization 4: Final tuning to make sure the coil-scalp distance is
        % equal to the designated value
        f = @(x)coilpos_cost_4_closeup(x,matsimnibs_3,shp,coilpos_base,tmp_target_norm,skin_norm,coil_to_scalp_distance);
        options = optimset('MaxIter',100,'MaxFunEvals',100);
        zscale = 1;
        lb = -coil_to_scalp_distance*zscale;
        ub = coil_to_scalp_distance*zscale;
        [x,fval4] = fminbnd(f,lb,ub,options);
        while fval4>0.1 && zscale<=16
            zscale = zscale*4;
            lb0 = lb;
            ub0 = ub;
            lb = (lb-zscale)*(x<=lb0+(ub0-lb0)/100) + ub0*(x>=ub0-(ub0-lb0)/100);
            ub = (ub+zscale)*(x>=ub0-(ub0-lb0)/100) + lb0*(x<=lb0+(ub0-lb0)/100);
            [x,fval4] = fminbnd(f,-coil_to_scalp_distance*zscale*(x<0),coil_to_scalp_distance*zscale*(x>0),options);
        end
		
		zshift_4 = x(1);
		if fval4<fval3
			coil_center_4 = coil_center_3+tmp_target_norm*zshift_4; % Ensure no intrusion at first
		else
			coil_center_4 = coil_center_3;
			fval4 = fval3;
		end
		matsimnibs_final = matsimnibs_3;
		matsimnibs_final(1:3,4) = coil_center_4(:);
		D_ROI = norm(target_ROI - coil_center_4);
		fprintf('Cost = %d -> %d -> %d -> %d\n',fval1,fval2,fval3,fval4);
    else
        matsimnibs_final = matsimnibs_2;
		matsimnibs_final(1:3,4) = coil_center_2(:);
        D_ROI = norm(target_ROI - coil_center_2);
        fprintf('Cost = %d -> %d\n',fval1,fval2);
    end
    
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
    D_skin = min(minD_valid(minD_valid>0));
    if min(abs(minD_true(angle_list<0)))<coil_to_scalp_distance*2
        inShapeCost = sum(angle_list<0);
    else
        inShapeCost = 0;
    end

    if isempty(D_skin)
        D_skin = min(minD_list);
    end

    flipCoil = false;
    tmpCR = target_ROI'-matsimnibs_final(1:3,4);
    tmpCR = tmpCR./norm(tmpCR);
    tmpCD = matsimnibs_final(1:3,1:3)*[0;0;1];
    if atan2(norm(cross(tmpCR,tmpCD)), dot(tmpCR,tmpCD))>pi/6
        flipCoil = true;
    end
    
    if inShapeCost>0 || flipCoil || abs(D_skin-coil_to_scalp_distance)>0.1
		matsimnibs_final = zeros(4,4);
		D_ROI = Inf;
		D_skin = Inf;
    end
end
return
