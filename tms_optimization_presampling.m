%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Generate spatially constrained TMS coil candidate placements for SimNIBS optimization     %
% Dependencies: SimNIBS-4.1, SPM12 (for custom defined ROI)                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;

% Determine which simulation(s) to run
protocol_pre_num = 1;
protocol_post_num = 1;
coil_num = 1;
brain_num = 1;

protocol_pre_list = {'direct','extension','direct_FT','extension_FT'};
protocol_post_list = {'cortex','cerebellum','inion'};
coil_types = {'D70','120BFV','DB80','DCC'}; % not included in the current study
coil_files = {'MagStim_D70.stl','Deymed_120BFV.stl','MagVenture_Cool-D-B80.stl','MagStim_DCC.stl'};

% SIMNIBSDIR = '~/SimNIBS-4.1/simnibs_env/lib/python3.9/site-packages/simnibs';
% addpath('~/matlab/spm12'); % Change to actual path to spm12

addpath(fullfile(SIMNIBSDIR,'matlab_tools'));
addpath(pwd);

cur_dir = pwd;
atlas_dir = pwd;
par_on = 0; % Enable parpool = 1, else = 0

%% Parameters
fine_tuning = false; % Fine tuning of the rotation matrix (slightly slower)
fine_collision_detection = true; % false: vertex-based detection; true: face-based detection (more robust but computationally expensive)

min_num_of_coil_positions = 60; % Mimimum number of points to evaluate around ROI
max_num_of_coil_positions = 80; % Maximum number of points to evaluate around ROI
coil_to_scalp_distance = 4; % Scalp-coil distance (in mm)
search_radius = 5; % Radius (in mm) of the sphere within which coil positions are searched
spatial_resolution = 1; % Average spatial resolution (in mm) of coil positions to be searched
search_angle = 360; % Angle range (in degree) of coil orientations to be searched
angle_resolution = 10; % Resolution (in degree) of coil orientations to be searched

fixed_coilpos = false; % True if only a specific coil position is to be optimized (w.r.t. the orientation)
% The next two parameters must be defined if fixed_coilpos==true; otherwise they have no effect
relative_position_to_Iz = [0,0,0]; % e.g., 3L1I = [-30,0,-10]
define_coilpos_by_fiducial = true; % true if the coil position is defined by Iz; false if by a MNI coordinate

init_angle = 0; % Initial orientation of the coil (useful if only one or a specific range of angles are to be tested)
mesh_truncate_nonROI = false; % truncate head mesh far from the ROI for even faster computations
mesh_truncate_nonROI_radius = 100; % Threshold radius of the truncated non-ROI region (in addition to coil dimensions)
xangle_range = pi/4; % Range of coil orientations along the z axis to optimize from
yangle_range = pi/4; % Range of coil orientations along the z axis to optimize from
xyangle_tolerance = pi/180; % 1 deg tolerance
coil_depth_threshold = 0.8; % Dismiss cases where coil-scalp distance is larger than coil_depth*coil_depth_threshold

%% Set up parpool
if par_on==1
    delete(gcp('nocreate'));
    numCPUs = feature('numCores');
    fprintf('Number of CPUs requested = %g\n',numCPUs);
    if numCPUs > 4 % on cluster
        pc_storage_dir = fullfile('pc_storage',getenv('SLURM_JOB_ID'));
        mkdir(pc_storage_dir);
        pc = parcluster('local');
        pc.JobStorageLocation =  pc_storage_dir;
    else
        pc = parcluster('local'); % use default
    end
    poolobj = parpool(pc,numCPUs);

    % Add cost function files
    addAttachedFiles(poolobj,'coilpos_optimization.m');
    addAttachedFiles(poolobj,'coilpos_cost_1_xyangle.m');
    addAttachedFiles(poolobj,'coilpos_cost_2_zshift.m');
    addAttachedFiles(poolobj,'coilpos_cost_3_finetune.m');
    addAttachedFiles(poolobj,'coilpos_cost_4_closeup.m');
    addAttachedFiles(poolobj,'simplifiedIntersect.m');
    addAttachedFiles(poolobj,'check_collision.m');
end

protocol_name = strcat(protocol_pre_list{protocol_pre_num},'_',protocol_post_list{protocol_post_num});

% Target ROI coordinates in MNI corresponding to both protocols
if contains(protocol_name,'cortex')
    % target_ROI_MNI = [-37, -21, 58];
    target_ROI_all = [...
    -5.4973   44.0089  214.4483
    0.6973   48.3736  197.7475
    -6.1587   32.5336  190.0981
    -7.0679   37.3732  167.1733
    -2.8196   54.6488  224.4266...
    ];
elseif contains(protocol_name,'cerebellum')
    % target_ROI_MNI = [-21, -56, -54];
    target_ROI_all = [...
    10.4586   15.6250  116.1754
    8.5062   10.1185  101.7312
    8.2875   27.0561   87.6023
    16.0376   13.9199   73.3336
    8.2679   22.0993  126.4091...
    ];
elseif contains(protocol_name,'inion')
    % target_ROI_MNI = [0, -46, -48];
    target_ROI_all = [...
    29.0211   24.1447  123.8895
    26.8520   16.7249  105.2721
    28.1902   34.1921   96.4108
    33.7427   23.6557   80.6734
    27.6133   29.2660  130.6377
    ];
end

% Subject folder
cd(strcat(atlas_dir,filesep,'brain',num2str(brain_num),'_charm'));
subpath = strcat('m2m_brain',num2str(brain_num),'_charm');

datadir_atlas = pwd;

% Load T1 image info
fname_t1 = ['brain' num2str(brain_num) '_t1_extended.nii'];
Data2Read_t1 = fullfile(datadir_atlas,fname_t1);
HeaderInfo_t1 = spm_vol(Data2Read_t1);
% Only for standard Sform (upper triangular matrix); otherwise change the following accordingly
xorigin = -HeaderInfo_t1.mat(1,4)/HeaderInfo_t1.mat(1,1);
yorigin = -HeaderInfo_t1.mat(2,4)/HeaderInfo_t1.mat(2,2);
zorigin = -HeaderInfo_t1.mat(3,4)/HeaderInfo_t1.mat(3,3);
coord_origin = round([xorigin,yorigin,zorigin]);
pixdim = [HeaderInfo_t1.mat(1,1) HeaderInfo_t1.mat(2,2) HeaderInfo_t1.mat(3,3)];

%% Presample feasible coil placements
%% Load head mesh
% Read the simulation result
head_mesh = mesh_load_gmsh4(...
    strcat(subpath,filesep,'brain',num2str(brain_num),'_charm','.msh') ...
    );

% Keep only skin surface (tag [5, 1005] in the mesh)
% skin_msh = mesh_extract_regions(head_mesh, 'region_idx', [5, 1005]);
skin_msh = mesh_extract_regions(head_mesh, 'region_idx', 1005);

% Load segmentation image
Data2Read_seg = fullfile(datadir_atlas,strcat('m2m_brain',num2str(brain_num),'_charm'),'final_tissues.nii.gz');
HeaderInfo_seg = spm_vol(Data2Read_seg);
[Data_seg,~] = spm_read_vols(HeaderInfo_seg);
Data_seg = Data_seg>0;

% Fill the head (vectorized where possible)
for iii = 1:size(Data_seg,3)
    for jjj = 1:size(Data_seg,1)
        rangeind = find(Data_seg(jjj,:,iii)>0);
        if length(rangeind)>2
            Data_seg(jjj,rangeind(1):rangeind(end),iii) = true;
        end
    end
end

se = strel('cube',9);
Data_seg_erode = imerode(Data_seg>0, se);

% Points constituting the head surface (vectorized shift calculation)
skin_msh_nodes = skin_msh.nodes;
skin_msh_nodes_shift = bsxfun(@plus, round(bsxfun(@rdivide, skin_msh_nodes, pixdim(1))), coord_origin);

% Generate exclude_pts in a more optimized way
exclude_pts = zeros(size(skin_msh_nodes_shift,1),1);

% Use linear indexing for better performance when possible
for iii = 1:length(exclude_pts)
    x = skin_msh_nodes_shift(iii,1);
    y = skin_msh_nodes_shift(iii,2);
    z = skin_msh_nodes_shift(iii,3);
    
    % Check bounds to avoid index errors
    if x >= 1 && x <= size(Data_seg_erode,1) && ...
       y >= 1 && y <= size(Data_seg_erode,2) && ...
       z >= 1 && z <= size(Data_seg_erode,3)
        exclude_pts(iii) = Data_seg_erode(x,y,z);
    end
end

include_pts = setdiff(1:size(skin_msh_nodes,1),find(exclude_pts>0));

% Update the triangles to use the filtered nodes
% Create a mapping from old indices to new indices
old_to_new = zeros(size(skin_msh.nodes, 1), 1);
old_to_new(include_pts) = 1:length(include_pts);

% Filter the triangles to keep only those with all vertices in include_pts (vectorized)
skin_triangles = skin_msh.triangles;
v1_valid = ismember(skin_triangles(:,1), include_pts);
v2_valid = ismember(skin_triangles(:,2), include_pts);
v3_valid = ismember(skin_triangles(:,3), include_pts);
valid_triangles = v1_valid & v2_valid & v3_valid;
skin_triangles = skin_triangles(valid_triangles, :);

% Update the triangle indices (vectorized)
skin_triangles(:,1) = old_to_new(skin_triangles(:,1));
skin_triangles(:,2) = old_to_new(skin_triangles(:,2));
skin_triangles(:,3) = old_to_new(skin_triangles(:,3));

% Update nodes
skin_msh_nodes = skin_msh.nodes(include_pts,:);

% Create a new skin mesh structure
skin_msh_filtered = struct();
skin_msh_filtered.nodes = skin_msh_nodes;
skin_msh_filtered.triangles = skin_triangles;

% Create alphaShape for the filtered skin mesh
shp = alphaShape(skin_msh_nodes, Inf);

clear mesh_downsample_index head_mesh Data2Read_seg HeaderInfo_seg Data_seg Data_seg_erode

%% Calculate surface normal vectors of the ROI
skin_norm = triangle_normals(skin_triangles, skin_msh_nodes, 20);
skin_norm(isnan(skin_norm))=0;

% Calculate head center once (vectorized min/max)
head_center = [mean([max(skin_msh_nodes(:,1)),min(skin_msh_nodes(:,1))]),...
    mean([max(skin_msh_nodes(:,2)),min(skin_msh_nodes(:,2))]),...
    mean([max(skin_msh_nodes(:,3)),min(skin_msh_nodes(:,3))])];

% Vectorized version of normal direction calculation
tmpvec1 = sqrt(sum((bsxfun(@minus, head_center, (skin_msh_nodes + skin_norm))).^2, 2));
tmpvec2 = sqrt(sum((bsxfun(@minus, head_center, (skin_msh_nodes - skin_norm))).^2, 2));
flip_mask = tmpvec1-tmpvec2 < 0;
skin_norm(flip_mask,:) = -skin_norm(flip_mask,:);

% Recalculate normals at neck bottom
neckbottom_index = find(min(skin_msh_nodes(:,3))-skin_msh_nodes(:,3)<=40);

if ~isempty(neckbottom_index)
    N_neckbottom = skin_norm(neckbottom_index,:);
    
    % Create reference point for neck normals
    neck_ref_point = head_center + [0, 0, min(skin_msh_nodes(:,3))-head_center(3)+50];
    
    % Vectorized neck normal calculation
    tmpvec1 = sqrt(sum((bsxfun(@minus, neck_ref_point, (skin_msh_nodes(neckbottom_index,:) + N_neckbottom))).^2, 2));
    tmpvec2 = sqrt(sum((bsxfun(@minus, neck_ref_point, (skin_msh_nodes(neckbottom_index,:) - N_neckbottom))).^2, 2));
    flip_neck_mask = tmpvec1-tmpvec2 < 0;
    N_neckbottom(flip_neck_mask,:) = -N_neckbottom(flip_neck_mask,:);
    skin_norm(neckbottom_index,:) = N_neckbottom;
end

clear exclude_pts N_neckbottom tmpvec1 tmpvec2 neckbottom_index

skin_norm_tmp = skin_norm;

%% Define the ROI
if ~isempty(protocol_name)
    protocol_name_sub = strcat('_',protocol_name);
else
    protocol_name_sub = protocol_name;
end

% ROI is a given coordinate
target_ROI = target_ROI_all(brain_num,:);

custom_target_coords_file = fullfile(pwd,strcat('brain',num2str(brain_num),'_target_label',protocol_name_sub,'.csv'));
dlmwrite(custom_target_coords_file,target_ROI);

% Load prespecified coil position
if fixed_coilpos % If only optimize a specific coil position w.r.t. the orientation
    if define_coilpos_by_fiducial
        % Generate a coil position based on the location of fiducials
        tmpf = fullfile(pwd,strcat('m2m_brain',num2str(brain_num),'_charm'),'eeg_positions','Fiducials.csv');
        fileID = fopen(tmpf,'r');
        fgetl(fileID);
        Iz_file = fgetl(fileID);
        tmpstr = strfind(Iz_file,',');
        xshift = str2double(Iz_file(tmpstr(1)+1:tmpstr(2)-1));
        yshift = str2double(Iz_file(tmpstr(2)+1:tmpstr(3)-1));
        zshift = str2double(Iz_file(tmpstr(3)+1:tmpstr(4)-1));
        fclose(fileID);
        Iz = [xshift,yshift,zshift];
        coilpos_init = Iz + relative_position_to_Iz;
    else
        % Pick an MNI coordinate
        coilpos_init = mni2subject_coords([-37, -21, 58], strcat('m2m_brain',num2str(brain_num),'_charm'));
    end
    
    % Find the closest point on the skin mesh to the coil position using alphaShape
    [closest_idx, ~] = nearestNeighbor(shp, coilpos_init);
    targetpos_skin = skin_msh_nodes(closest_idx,:);
    target_norm = skin_norm_tmp(closest_idx,:);
else
    % Locate projections of the ROI on head skin based on distance
    search_radius_tmp = search_radius;
    
    % Vectorized distance calculation
    distances = sqrt(sum(bsxfun(@minus, skin_msh_nodes, target_ROI).^2, 2));
    targetpos_skin_ind = find(distances <= search_radius_tmp);
    targetpos_skin = skin_msh_nodes(targetpos_skin_ind,:);

    % Ensure we have enough coil positions
    while isempty(targetpos_skin_ind) || length(targetpos_skin_ind) < min_num_of_coil_positions
        search_radius_tmp = search_radius_tmp + spatial_resolution;
        targetpos_skin_ind = find(distances <= search_radius_tmp);
        targetpos_skin = skin_msh_nodes(targetpos_skin_ind,:);
        
        if ~isempty(targetpos_skin_ind) && length(targetpos_skin_ind) >= min_num_of_coil_positions
            fprintf('Searching radius too small. Increased to %d mm.\n', search_radius_tmp);
        end
    end

    % Sort by distance from the ROI center (vectorized)
    [~, idx_sort] = sort(distances(targetpos_skin_ind));
    max_positions = min([max_num_of_coil_positions, length(idx_sort)]);
    targetpos_skin = targetpos_skin(idx_sort(1:max_positions),:);
    target_norm = skin_norm_tmp(targetpos_skin_ind(idx_sort(1:max_positions)),:);

    % Handle zero-value normal vectors (vectorized)
    zeronormind = sum(target_norm==0,2)==3;
    validnormind = sum(target_norm==0,2)~=3;
    if sum(validnormind) > 0
        target_norm_mean = mean(target_norm(validnormind,:), 1);
        if norm(target_norm_mean) > 0
            target_norm(zeronormind,:) = repmat(target_norm_mean/norm(target_norm_mean), sum(zeronormind), 1);
        end
    end
end

% Truncate the skin mesh if requested (to reduce computation time)
if mesh_truncate_nonROI
    if size(targetpos_skin,1)>1
        tmpSkinTarget = mean(targetpos_skin, 1); % Vectorized mean
    else
        tmpSkinTarget = targetpos_skin;
    end
    tmpOppositeTarget = head_center*2-tmpSkinTarget;
    
    % Find vertices to keep (exclude those far from target) - vectorized
    distances_to_opposite = sqrt(sum(bsxfun(@minus, skin_msh_nodes, tmpOppositeTarget).^2, 2));
    far_vertices = distances_to_opposite > mesh_truncate_nonROI_radius;
    
    % Create a new vertices list with only the included vertices
    skin_msh_nodes_filtered = skin_msh_nodes(far_vertices,:);
    skin_norm_filtered = skin_norm_tmp(far_vertices,:);
    
    % Update the triangles for the filtered mesh (vectorized)
    old_to_new = zeros(size(skin_msh_nodes, 1), 1);
    old_to_new(far_vertices) = 1:sum(far_vertices);
    
    % Filter triangles vectorized
    far_vertices_idx = find(far_vertices);
    v1_far = ismember(skin_triangles(:,1), far_vertices_idx);
    v2_far = ismember(skin_triangles(:,2), far_vertices_idx);
    v3_far = ismember(skin_triangles(:,3), far_vertices_idx);
    valid_triangles = v1_far & v2_far & v3_far;
    skin_triangles_filtered = skin_triangles(valid_triangles, :);
    
    % Update the triangle indices (vectorized)
    skin_triangles_filtered(:,1) = old_to_new(skin_triangles_filtered(:,1));
    skin_triangles_filtered(:,2) = old_to_new(skin_triangles_filtered(:,2));
    skin_triangles_filtered(:,3) = old_to_new(skin_triangles_filtered(:,3));
    
    % Update the mesh structure
    skin_msh_filtered.nodes = skin_msh_nodes_filtered;
    skin_msh_filtered.triangles = skin_triangles_filtered;
    
    skin_msh_nodes = skin_msh_nodes_filtered;
    skin_triangles = skin_triangles_filtered;
    skin_norm_tmp = skin_norm_filtered;
end

% Save mesh data for this brain
save(strcat('brain',num2str(brain_num),protocol_name_sub,'_msh.mat'),'skin_msh_filtered','target_norm','targetpos_skin','skin_norm');

clear zeronormind validnormind target_norm_mean idx_sort tmpSkinTarget tmpOppositeTarget

%% Iterate through each coil
% Select the coil model
coil_type = coil_types{coil_num};
coil_file = coil_files{coil_num};
fnamecoil = fullfile(SIMNIBSDIR,'resources','coil_models','Drakaki_BrainStim_2022',coil_file);
custom_coil_matrices_file = fullfile(pwd,strcat('matrices_brain',num2str(brain_num),protocol_name_sub,'_',coil_type,'_all.csv'));

% Read coil geometry - now including coilpos_face
load(fullfile(cur_dir,strcat('coilpos_base_',coil_type,'.mat')),'coilpos_base','coilpos_face');

%% Iterate through each target point
% Parameters
minx0 = min(abs(coilpos_base(:,1)));
Ix0 = abs(abs(coilpos_base(:,1))-minx0)<1;
maxz0 = max(coilpos_base(Ix0,3));
maxzl = max(coilpos_base(coilpos_base(:,1)>minx0,3));
maxzr = max(coilpos_base(coilpos_base(:,1)<-minx0,3));
% rangex = range(coilpos_base(:,1));
zrange = round((mean([maxzl,maxzr])-maxz0)*coil_depth_threshold+maxz0); % Coil depth

% Use pre-calculated thetaz values and adjust for D70 if needed
thetaz_all_base = linspace(-search_angle/2+init_angle, search_angle/2+init_angle, search_angle/angle_resolution+1)*pi/180;
thetaz_all = thetaz_all_base(1:end-1);
if strcmp(coil_type,'D70')
    thetaz_all = thetaz_all + pi; % D70 current facing up by default
end

thetaz_size = length(thetaz_all);

% Pre-allocate matrices for results
matrices_all = cell(thetaz_size, size(targetpos_skin,1));
D_ROI_all = zeros(thetaz_size, size(targetpos_skin,1));
D_skin_all = zeros(thetaz_size, size(targetpos_skin,1));

fprintf('Evaluating %d coil placements in total...\n', numel(D_ROI_all));

%% Outer loop - iterate through each skin vertex (which the coil center faces) above the ROI
tic;
if par_on==0
    for target_skin_idx = 1:size(targetpos_skin,1)
	    target_skin = targetpos_skin(target_skin_idx,:);
	    tmp_target_norm = target_norm(target_skin_idx,:);

        % Function to calculate rotation matrix from [0,0,-1] (original coil -zdir)
        if isequal(tmp_target_norm, [0,0,-1])
            U = eye(3,3);
        elseif isequal(tmp_target_norm, [0,0,1])
            U = [1,0,0;0,-1,0;0,0,-1];
        else
            ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
            RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
            U = RU([0,0,-1], tmp_target_norm);
        end

        %% Inner loop - rotate the coil along its own z-axis
        % Find the best x/y-axis rotations and shift along the normal vector
        for thetaz_idx = 1:thetaz_size
            thetaz = thetaz_all(thetaz_idx);
            if ~strcmp(coil_type,'D70') && abs(thetaz)>=pi/4
                fprintf('Cost = Inf\n');
                matsimnibs_final = zeros(4,4);
                D_ROI = Inf;
                D_skin = Inf; % D70 current facing up by default
            else
                % Pass coilpos_face to the optimization function
                if ~fine_collision_detection
                    [matsimnibs_final,D_ROI,D_skin] = coilpos_optimization(shp,coilpos_base,U,thetaz,maxz0,zrange,...
					    target_ROI,target_skin,skin_norm,tmp_target_norm,coil_to_scalp_distance,xangle_range,yangle_range,xyangle_tolerance,fine_tuning);
                else
                    [matsimnibs_final,D_ROI,D_skin] = coilpos_optimization_face(skin_msh_filtered, coilpos_base, coilpos_face, U, thetaz, maxz0, zrange,...
                        target_ROI, target_skin, skin_norm_tmp, tmp_target_norm, coil_to_scalp_distance, xangle_range, yangle_range, xyangle_tolerance, fine_tuning);
                end
            end
            matrices_all{thetaz_idx,target_skin_idx} = matsimnibs_final;
            D_ROI_all(thetaz_idx,target_skin_idx) = D_ROI;
            D_skin_all(thetaz_idx,target_skin_idx) = D_skin;
        end
        
        % Display progress periodically
        if mod(target_skin_idx, 10) == 0
            fprintf('Processed %d/%d target positions\n', target_skin_idx, size(targetpos_skin,1));
        end
    end
else
    % Parallel version
    % Pre-allocate storage for results from parallel workers
    local_matrices = cell(thetaz_size, size(targetpos_skin,1));
    local_D_ROI = zeros(thetaz_size, size(targetpos_skin,1));
    local_D_skin = zeros(thetaz_size, size(targetpos_skin,1));
    
    parfor target_skin_idx = 1:size(targetpos_skin,1)
	    target_skin = targetpos_skin(target_skin_idx,:);
	    tmp_target_norm = target_norm(target_skin_idx,:);

        % Function to calculate rotation matrix from [0,0,-1] (original coil -zdir)
        if isequal(tmp_target_norm, [0,0,-1])
            U = eye(3,3);
        elseif isequal(tmp_target_norm, [0,0,1])
            U = [1,0,0;0,-1,0;0,0,-1];
        else
            ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
            RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
            U = RU([0,0,-1], tmp_target_norm);
        end
        
        % Process all angles for this target position
        target_matrices = cell(thetaz_size, 1);
        target_D_ROI = zeros(thetaz_size, 1);
        target_D_skin = zeros(thetaz_size, 1);

        % Inner loop - rotate the coil along its own z-axis
        for thetaz_idx = 1:thetaz_size
            thetaz = thetaz_all(thetaz_idx);
            if ~strcmp(coil_type,'D70') && abs(thetaz)>=pi/4
                if thetaz_idx == 1 % Reduce console output
                    fprintf('Cost = Inf\n');
                end
                matsimnibs_final = zeros(4,4);
                D_ROI = Inf;
                D_skin = Inf; % D70 current facing up by default
            else
                % Pass coilpos_face to the optimization function
                if ~fine_collision_detection
                    [matsimnibs_final,D_ROI,D_skin] = coilpos_optimization(shp,coilpos_base,U,thetaz,maxz0,zrange,...
					    target_ROI,target_skin,skin_norm,tmp_target_norm,coil_to_scalp_distance,xangle_range,yangle_range,xyangle_tolerance,fine_tuning);
                else
                    [matsimnibs_final,D_ROI,D_skin] = coilpos_optimization_face(skin_msh_filtered, coilpos_base, coilpos_face, U, thetaz, maxz0, zrange,...
                        target_ROI, target_skin, skin_norm_tmp, tmp_target_norm, coil_to_scalp_distance, xangle_range, yangle_range, xyangle_tolerance, fine_tuning);
                end
            end
            target_matrices{thetaz_idx} = matsimnibs_final;
            target_D_ROI(thetaz_idx) = D_ROI;
            target_D_skin(thetaz_idx) = D_skin;
        end
        
        % Store results for this target position
        for thetaz_idx = 1:thetaz_size
            local_matrices{thetaz_idx, target_skin_idx} = target_matrices{thetaz_idx};
            local_D_ROI(thetaz_idx, target_skin_idx) = target_D_ROI(thetaz_idx);
            local_D_skin(thetaz_idx, target_skin_idx) = target_D_skin(thetaz_idx);
        end
    end
    
    % Copy results back to main variables
    matrices_all = local_matrices;
    D_ROI_all = local_D_ROI;
    D_skin_all = local_D_skin;
end
toc;

% Save matrices files
dotind = strfind(custom_coil_matrices_file,'.');
save(strcat(custom_coil_matrices_file(1:dotind(end)),'mat'),'matrices_all','D_ROI_all','D_skin_all');

cd(cur_dir);
if par_on==1
    delete(poolobj);
    if numCPUs > 4
        rmdir(pc_storage_dir,'s');
    end
end

fprintf('TMS coil optimization completed\n');
