%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Implement SimNIBS optimization of spatially constrained candidate TMS coil placements     %
% Compatible with both Direct (w/ or w/o pardiso solver) and ADM                            %
% Dependencies: SimNIBS-4.1                                                                 %
% Requires modified SimNIBS scripts (by Xu Zhang) to work: opt_struct.py, optimize_tms.py   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;

protocol_pre_num = 1;
protocol_post_num = 1;
brain_num = 1;
coil_num = 1;

% SIMNIBSDIR = '~/SimNIBS-4.1/simnibs_env/lib/python3.9/site-packages/simnibs';
% addpath('~/matlab/spm12'); % Change to actual path to spm12

addpath(fullfile(SIMNIBSDIR,'matlab_tools'));
addpath(fullfile(pwd));

% Specify the optimization method (ADM or direct)
tms_opt = opt_struct('TMSoptimize');

%%% Run ADM:
tms_opt.method = 'ADM';
tms_opt.resampleflag = false; % Disable ADM resampling of coil geometry

%%% Run direct method:
% tms_opt.method = 'direct';
% tms_opt.solver_options = 'pardiso';

% unsetenv('LD_LIBRARY_PATH'); % SimNIBS 4.1 may not run properly without this command in Unix

protocol_pre_list = {'direct','extension','direct_FT','extension_FT'};
protocol_post_list = {'cortex','cerebellum','inion'};
coil_types = {'D70','120BFV','DB80','DCC'}; % not included in the current study
coil_files = {'MagStim_D70.stl','Deymed_120BFV.stl','MagVenture_Cool-D-B80.stl','MagStim_DCC.stl'};

cur_dir = fullfile(pwd);
atlas_dir = fullfile(pwd);

protocol_name = strcat(protocol_pre_list{protocol_pre_num},'_',protocol_post_list{protocol_post_num});

disp(brain_num);
cd(strcat(atlas_dir,filesep,'brain',num2str(brain_num),'_charm'));

coil_type = coil_types{coil_num};
coil_file = coil_files{coil_num};

if ~isempty(tms_opt.solver_options)
	solver_options = strcat('_',tms_opt.solver_options);
else
	solver_options = '';
end

tms_opt.open_in_gmsh = false;

% Subject folder
tms_opt.subpath = fullfile(atlas_dir,strcat('brain',num2str(brain_num),'_charm'),strcat('m2m_brain',num2str(brain_num),'_charm'));
% Select output folder
tms_opt.pathfem = fullfile(atlas_dir,strcat('brain',num2str(brain_num),'_',protocol_name,'_ADM_',coil_type));

if exist(strcat(tms_opt.pathfem,'_simulation_interp'),'dir')
	disp('Already done.');
	quit();
end

% Scalp-coil distance (minimum distance if custom_coil_matrices_option == true)
tms_opt.distance = 4;

% Bypassed if tms_opt.custom_coil_matrices_option == true;
tms_opt.search_radius = 20;
tms_opt.spatial_resolution = 2;
tms_opt.search_angle = 0;
tms_opt.angle_resolution = 10;

% Non-SimNIBS parameters
ROI_spatial_resolution = 5;

tms_opt.fnamecoil = fullfile(SIMNIBSDIR,'resources','coil_models','Drakaki_BrainStim_2022',coil_file);

% Read coil geometry
load(fullfile(cur_dir,strcat('coilpos_base_',coil_type,'.mat')),'coilpos_base');

tms_opt.custom_target_coords_option = false;
tms_opt.custom_coil_matrices_option = true;
tms_opt.custom_target_coords_file = fullfile(pwd,strcat('brain',num2str(brain_num),'_target_label_',protocol_name,'.csv'));
tms_opt.custom_coil_matrices_file = fullfile(pwd,strcat('matrices_brain',num2str(brain_num),'_',protocol_name,'_',coil_type,'_all.csv'));

% Readjust spatial resolution of ROI
ROI_subject = dlmread(tms_opt.custom_target_coords_file);
if tms_opt.custom_target_coords_option % Target coordinates are fully provided
	if size(ROI_subject,1)>1 % If need to downsample
		ROI_subject = unique(round(ROI_subject/ROI_spatial_resolution)*ROI_spatial_resolution,'rows');
	end
	if size(ROI_subject,1)==1 % Duplicate rows so the subsequent Python scheme can take the average over rows
		ROI_subject = repmat(ROI_subject,4,1);
	end
	tms_opt.custom_target_coords_file = fullfile(pwd,strcat('brain',num2str(brain_num),'_target_label_downsampled_',protocol_name,'.csv'));
	dlmwrite(tms_opt.custom_target_coords_file,ROI_subject);
	tms_opt.target = mean(dlmread(tms_opt.custom_target_coords_file));
else % Only for testing; the target coordinate(s) should be specified in tms_optimization_presampling.m
	% Target Cerebellar lobule VIII (L); note that SimNIBS may project this ROI onto the ear canal
	tms_opt.target = ROI_subject;
	tms_opt.target_size = 5;
end

dotind = strfind(tms_opt.custom_coil_matrices_file,'.');
custom_coil_matrices = load(strcat(tms_opt.custom_coil_matrices_file(1:dotind(end)),'mat'));

matrices_all = reshape(custom_coil_matrices.matrices_all,size(custom_coil_matrices.matrices_all,1)*size(custom_coil_matrices.matrices_all,2),1);
D_skin_all = reshape(custom_coil_matrices.D_skin_all,size(custom_coil_matrices.D_skin_all,1)*size(custom_coil_matrices.D_skin_all,2),1);
D_ROI_all = reshape(custom_coil_matrices.D_ROI_all,size(custom_coil_matrices.D_ROI_all,1)*size(custom_coil_matrices.D_ROI_all,2),1);

% Only run feasible coil placements
minx0 = min(abs(coilpos_base(:,1)));
Ix0 = find(abs(coilpos_base(:,1))==minx0);
maxz0 = max(coilpos_base(Ix0,3));
maxzl = max(coilpos_base(coilpos_base(:,1)>minx0,3));
maxzr = max(coilpos_base(coilpos_base(:,1)<minx0,3));
zrange = round((max(maxzl,maxzr)-maxz0)*0.4+maxz0); % Coil depth
If_exclude = find(abs(D_skin_all-tms_opt.distance)>0.1);
matrices_all(If_exclude) = [];
D_skin_all(If_exclude) = [];
D_ROI_all(If_exclude) = [];

ind_all = length(D_ROI_all(:));
mat_all = zeros(ind_all*4,4);
for i = 1:ind_all
	mat_all((i*4-3):(i*4),:) = matrices_all{i};
	mat_all((i*4),4) = 1;
end
dlmwrite(tms_opt.custom_coil_matrices_file,mat_all,'precision','%.15f');

% Run optimization
try
	run_simnibs(tms_opt);
catch
	quit();
end

% Re-run simulation with volume options
tmpf = dir(fullfile(strcat(tms_opt.pathfem)));
for tff = 3:length(tmpf)
	if strcmp(tmpf(tff).name(end-2:end),'log')
		optresult = fullfile(tms_opt.pathfem,tmpf(tff).name);
		break;
	end
end

fileID = fopen(optresult,'r');
tmpline = fgetl(fileID);
while ~strcmp(tmpline(1),'=')
	tmpline = fgetl(fileID);
end
matline1 = fgetl(fileID);
matline2 = fgetl(fileID);
matline3 = fgetl(fileID);
matline4 = fgetl(fileID);
matline_all = [matline1(2:end-1) ';' matline2(3:end-1) ';' matline3(3:end-1) ';' matline4(3:end-1) ';' ];
matsimnibs = eval(matline_all);
fclose(fileID);

S = sim_struct('SESSION');
%% Run SIMNIBS
S.subpath = tms_opt.subpath;
S.poslist{1} = sim_struct('TMSLIST');
S.poslist{1}.fnamecoil = tms_opt.fnamecoil;
S.poslist{1}.pos(1).matsimnibs = matsimnibs;

% Simulation with interpolation
% disp('Re-running simulation with interpolated volumes...');
S.pathfem = strcat(tms_opt.pathfem,'_simulation_interp');
S.open_in_gmsh=false; % Open result in gmsh when done
S.map_to_vol=false; % write fields as nifti
S.map_to_MNI=false; % write fields as nifti in MNI space
S.tissues_in_niftis = 'all';  % determines for which tissues the fields will be written
S.fiducials = sim_struct('FIDUCIALS');
run_simnibs(S);
cd(cur_dir);
