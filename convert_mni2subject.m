% Example script to run a TMS optimization
% Guilherme Saturnino, 2019

clear;clc;
% addpath('~/SimNIBS-4.0/simnibs_env/lib/python3.9/site-packages/simnibs/matlab_tools');
% tms_opt = opt_struct('TMSoptimize');

for brain_num = 1:5
    % Subject folder
    subpath = strcat('brain',num2str(brain_num),'_charm\m2m_brain',num2str(brain_num),'_charm');
    % Select a target for the optimization
%     disp(mni2subject_coords([-37, -21, 58], subpath));
%     disp(mni2subject_coords([-21, -56, -54], subpath));
    disp(mni2subject_coords([0, -46, -48], subpath));
end

%% Cortex
% -5.4973   44.0089  214.4483
% 0.6973   48.3736  197.7475
% -6.1587   32.5336  190.0981
% -7.0679   37.3732  167.1733
% -2.8196   54.6488  224.4266

%% Cerebellum
% 10.4586   15.6250  116.1754
% 8.5062   10.1185  101.7312
% 8.2875   27.0561   87.6023
% 16.0376   13.9199   73.3336
% 8.2679   22.0993  126.4091

%% Inion
% 29.0211   24.1447  123.8895
% 26.8520   16.7249  105.2721
% 28.1902   34.1921   96.4108
% 33.7427   23.6557   80.6734
% 27.6133   29.2660  130.6377
