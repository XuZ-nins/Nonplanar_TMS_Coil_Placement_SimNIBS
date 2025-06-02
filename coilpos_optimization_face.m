function [matsimnibs_final,D_ROI,D_skin] = coilpos_optimization_face(skin_msh, coilpos_base, coilpos_face, U, thetaz, maxz0, zrange,...
    target_ROI, target_skin, ~, tmp_target_norm, coil_to_scalp_distance, xangle_range, yangle_range, xyangle_tolerance, fine_tuning)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Main function for coil position optimization                                              %
% Modified to use improved collision detection                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize transformation matrix
matsimnibs_orig = zeros(4,4);
matsimnibs_orig(4,4) = 1;

% Initialize coil far from the head to avoid initial collisions
coil_center_1 = target_skin+tmp_target_norm*(max(zrange+15,30)+coil_to_scalp_distance); 

% Create initial rotation matrix around z-axis
rotMatz_1 = [cos(thetaz),-sin(thetaz) 0;sin(thetaz),cos(thetaz),0;0,0,1];
matsimnibs_orig(1:3,1:3) = rotMatz_1;

% ---- Optimization Stage 1: Find best x/y-axis rotations for coil orientation ----
% Configure optimization options
options = optimoptions(@fmincon,'MaxFunctionEvaluations',1500,'StepTolerance',xyangle_tolerance, ...
    'Display','off', 'Algorithm','sqp');

% Create cost function for finding optimal x/y rotations
f = @(x)coilpos_cost_1_xyangle_face(x, matsimnibs_orig, coil_center_1, U, coilpos_base, coilpos_face, skin_msh, target_ROI, coil_to_scalp_distance);

% Set bounds and run optimization with multiple starting points
x0 = [xangle_range/6, yangle_range/6];
lb = [-xangle_range -yangle_range];
ub = [xangle_range yangle_range];

[x, fval1] = fmincon(f, x0, [], [], [], [], lb, ub, [], options);

% Use the optimal angles to create rotation matrices
thetax_1 = x(1); % Clockwise = positive
thetay_1 = x(2); % Clockwise = positive
rotMatx_1 = [1, 0, 0;
             0, cos(thetax_1), -sin(thetax_1);
             0, sin(thetax_1), cos(thetax_1)];
rotMaty_1 = [cos(thetay_1), 0, sin(thetay_1);
             0, 1, 0;
             -sin(thetay_1), 0, cos(thetay_1)];

% Update transformation matrix with rotations
matsimnibs_1 = matsimnibs_orig;
matsimnibs_1(1:3,1:3) = U*rotMatx_1*rotMaty_1*matsimnibs_1(1:3,1:3);
matsimnibs_1(1:3,4) = coil_center_1';
coil_vec_1 = transform_normal(matsimnibs_1, [0,0,1]);
ROI_vec_1 = target_ROI-coil_center_1; ROI_vec_1 = ROI_vec_1/norm(ROI_vec_1);

% Check for infeasible solution 
if fval1 > 0 || angle_between(ROI_vec_1,coil_vec_1)>30 % Coil focus deviates too much from ROI
    fprintf('Cost = Inf (Stage 1 failed)\n');
    matsimnibs_final = zeros(4,4);
    D_ROI = Inf;
    D_skin = Inf;
    return;
end

% ---- Optimization Stage 2: Find best shift along normal vector to minimize coil-scalp distance ----
% Create cost function for optimal z-shift
f = @(x)coilpos_cost_2_zshift_face(x, matsimnibs_1, coilpos_base, coilpos_face, skin_msh, target_skin, tmp_target_norm, coil_to_scalp_distance);

% Configure optimization options - optimized for faster convergence
options = optimset('MaxIter', 100, ...         % Reduced from 150
                   'MaxFunEvals', 100, ...     % Reduced from 150
                   'TolX', 0.2, ...            % Increased from 0.1 for faster convergence
                   'TolFun', 0.2, ...          % Added function tolerance
                   'Display', 'off');          % Disable output

% Set bounds for z-shift search - narrower range for faster convergence
lb = maxz0;
ub = max(zrange+15, 30) + maxz0;

% Run optimization
[x, fval2] = fminbnd(f, lb, ub, options);

% If optimization fails, try grid search with fewer points
if fval2 > 10
    % Use fewer points for faster grid search (50 instead of 100)
    num_points = 50;
    cost_all = zeros(num_points, 1);
    zshift_all = linspace(lb, ub, num_points);
    
    for ii = 1:num_points
        zshift_tmp = zshift_all(ii);
        cost_all(ii) = f(zshift_tmp);
    end
    
    [fval2, Iminc] = min(cost_all);
    x = zshift_all(Iminc);
end

% Store optimal z-shift and update transformation matrix
zshift_2 = x(1);
matsimnibs_2 = matsimnibs_orig;
matsimnibs_2(1:3,1:3) = U*rotMatx_1*rotMaty_1*matsimnibs_orig(1:3,1:3);
coil_center_2 = target_skin + tmp_target_norm*(zshift_2 + coil_to_scalp_distance);
matsimnibs_2(1:3,4) = coil_center_2(:);

% ---- Optimization Stage 3 (if enabled): Fine-tune x/y rotations and z shift together ----
if fine_tuning
    % Create cost function for simultaneous fine-tuning of angles and z-shift
    f = @(x)coilpos_cost_3_finetune_face(x, matsimnibs_orig, U, thetax_1, thetay_1, coil_center_2, target_skin, coilpos_base, coilpos_face, skin_msh, tmp_target_norm, coil_to_scalp_distance);
    
    % Configure optimization options with optimized parameters
    options = optimoptions(@fmincon,...
        'MaxFunctionEvaluations', 500,...        % Reduced from 2000
        'MaxIterations', 300,...                  % Added explicit iteration limit
        'StepTolerance', 1e-3,...
        'OptimalityTolerance', 1e-2,...           % Increased for faster convergence
        'FunctionTolerance', 1e-2,...             % Increased for faster convergence
        'Display', 'off',...                      % Changed from 'iter' to 'off'
        'Algorithm', 'sqp',...
        'ScaleProblem', 'obj-and-constr');        % Scale the problem

    % Set bounds for fine-tuning with better initial point
    x0 = [0, 0, 0];
    lb = [-xangle_range/3, -yangle_range/3, -coil_to_scalp_distance-zrange/4];
    ub = [xangle_range/3, yangle_range/3, coil_to_scalp_distance];
    
    [x, fval3] = fmincon(f, x0, [], [], [], [], lb, ub, [], options);
    
    % Extract and apply fine-tuned parameters
    thetax_3 = x(1); 
    thetay_3 = x(2); 
    zshift_3 = x(3);
    
    % Create fine-tuned rotation matrices
    rotMatx_3 = [1, 0, 0;
                 0, cos(thetax_1+thetax_3), -sin(thetax_1+thetax_3);
                 0, sin(thetax_1+thetax_3), cos(thetax_1+thetax_3)];
    rotMaty_3 = [cos(thetay_1+thetay_3), 0, sin(thetay_1+thetay_3);
                 0, 1, 0;
                 -sin(thetay_1+thetay_3), 0, cos(thetay_1+thetay_3)];
    
    % Update transformation matrix with fine-tuned parameters
    matsimnibs_3 = matsimnibs_orig;
    matsimnibs_3(1:3,1:3) = U*rotMatx_3*rotMaty_3*matsimnibs_orig(1:3,1:3);
    coil_center_3 = coil_center_2 + tmp_target_norm*zshift_3;
    matsimnibs_3(1:3,4) = coil_center_3(:);

    % ---- Optimization Stage 4: Final tuning to achieve exact coil-scalp distance ----
    % Create cost function focusing purely on achieving exact coil-scalp distance
    f = @(x)coilpos_cost_4_closeup_face(x, matsimnibs_3, coilpos_base, coilpos_face, skin_msh, tmp_target_norm, coil_to_scalp_distance);
    
    % Configure optimization options - optimized for faster convergence
    options = optimset('MaxIter', 80, ...          % Reduced from 150
                       'MaxFunEvals', 80, ...      % Reduced from 150
                       'TolX', 0.1, ...            % Increased from 0.05
                       'TolFun', 0.1, ...          % Added function tolerance
                       'Display', 'off');          % Disable output

    % Set initial scale and bounds - use narrower initial range
    zscale = 0.8;
    lb = -coil_to_scalp_distance*zscale;
    ub = coil_to_scalp_distance*zscale;
    
    % Run optimization
    [x, fval4] = fminbnd(f, lb, ub, options);
    
    % If optimization doesn't converge well, expand search range iteratively
    % Reduced max attempts and larger scaling for faster convergence
    attempt_count = 0;
    max_attempts = 3;  % Reduced from 5
    
    while fval4 > 0.15 && zscale <= 8 && attempt_count < max_attempts
        attempt_count = attempt_count + 1;
        zscale = zscale * 2;
        lb0 = lb;
        ub0 = ub;
        
        % Adjust bounds based on current solution
        if x <= lb0 + (ub0-lb0)/100
            lb = lb - zscale;
        elseif x >= ub0 - (ub0-lb0)/100
            ub = ub + zscale;
        else
            % If solution is in the middle, expand both bounds
            lb = lb - zscale/2;
            ub = ub + zscale/2;
        end
        
        % Run optimization with updated bounds
        [x, fval4] = fminbnd(f, lb, ub, options);
    end
    
    % Apply final z-shift adjustment
    zshift_4 = x(1);
    
    % Choose the best result from stages 3 and 4
    if fval4 < fval3
        coil_center_4 = coil_center_3 + tmp_target_norm*zshift_4;
    else
        coil_center_4 = coil_center_3;
        fval4 = fval3;
    end
    
    % Set final transformation matrix
    matsimnibs_final = matsimnibs_3;
    matsimnibs_final(1:3,4) = coil_center_4(:);
    D_ROI = norm(target_ROI - coil_center_4);
    fprintf('Cost = %d -> %d -> %d -> %d\n', fval1, fval2, fval3, fval4);
else
    % If no fine-tuning, use results from stage 2
    matsimnibs_final = matsimnibs_2;
    matsimnibs_final(1:3,4) = coil_center_2(:);
    D_ROI = norm(target_ROI - coil_center_2);
    fprintf('Cost = %d -> %d\n', fval1, fval2);
end

% Transform coil vertices using the final transformation matrix
coilpos = (matsimnibs_final*[coilpos_base(:,1:3), ones(size(coilpos_base,1),1)]')';
coilpos = coilpos(:,1:3);

% Get minimum distance between coil and skin for the final result
[D_skin, collision_flag] = get_coil_skin_distance(coilpos, coilpos_face, skin_msh, coil_to_scalp_distance);

if collision_flag>0
    matsimnibs_final = zeros(4,4);
    D_ROI = Inf;
    D_skin = 0;
end
