%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Enlarge the MRI space to make room for neck (and ear & nose) reconstruction               %
% Dependencies: SPM12                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all
% addpath('spm12'); % Change to actual path to spm12

for brain_num = 1:5

    rawdir = fullfile(pwd, 'raw');
    newdir = fullfile(pwd, strcat('brain',num2str(brain_num),'_charm'));
    
    %% Add neck to T1 images
    fname_brain = strcat('brain',num2str(brain_num),'_t1.nii');
    Data2Read_t1 = fullfile(rawdir,fname_brain);
    HeaderInfo_t1 = spm_vol(Data2Read_t1);
    Data_t1 = spm_read_vols(HeaderInfo_t1);
    
    % Downsample to 0.6 mm
    Data_t1_ds = imresize3(Data_t1,0.5,'linear');
    Data_t1_new = zeros(size(Data_t1_ds)+[100 100 250]);
    Data_t1_new(51:end-50,51:end-50,201:end-50) = Data_t1_ds;

    HeaderInfo_t1.fname = fullfile(newdir,replace(fname_brain,'.','_extended.'));
    HeaderInfo_t1.dim = size(Data_t1_new);
    HeaderInfo_t1.mat(1:3,1:3) = HeaderInfo_t1.mat(1:3,1:3)*2;
    spm_write_vol(HeaderInfo_t1,Data_t1_new);

    %% Add neck to T2 images
    fname_brain = strcat('brain',num2str(brain_num),'_t2.nii');
    Data2Read_t2 = fullfile(rawdir,fname_brain);
    HeaderInfo_t2 = spm_vol(Data2Read_t2);
    Data_t2 = spm_read_vols(HeaderInfo_t2);
    
    % Downsample to 0.6 mm
    Data_t2_ds = imresize3(Data_t2,0.5,'linear');
    Data_t2_new = zeros(size(Data_t2_ds)+[100 100 250]);
    Data_t2_new(51:end-50,51:end-50,201:end-50) = Data_t2_ds;

    HeaderInfo_t2.fname = fullfile(newdir,replace(fname_brain,'.','_extended.'));
    HeaderInfo_t2.dim = size(Data_t2_new);
    HeaderInfo_t2.mat(1:3,1:3) = HeaderInfo_t2.mat(1:3,1:3)*2;
    spm_write_vol(HeaderInfo_t2,Data_t2_new);

    %% Add neck to mask images
    fname_brain = strcat('brain',num2str(brain_num),'_labels.nii');
    Data2Read_labels = fullfile(rawdir,fname_brain);
    HeaderInfo_labels = spm_vol(Data2Read_labels);
    Data_labels = spm_read_vols(HeaderInfo_labels);
    
    % Downsample to 0.6 mm
    Data_labels_ds = imresize3(Data_labels,0.5,'nearest');
    Data_labels_new = zeros(size(Data_labels_ds)+[100 100 250]);
    Data_labels_new(51:end-50,51:end-50,201:end-50) = Data_labels_ds;

    HeaderInfo_labels.fname = fullfile(newdir,replace(fname_brain,'.','_extended.'));
    HeaderInfo_labels.dim = size(Data_labels_new);
    HeaderInfo_labels.mat(1:3,1:3) = HeaderInfo_labels.mat(1:3,1:3)*2;
    spm_write_vol(HeaderInfo_labels,Data_labels_new);
end
