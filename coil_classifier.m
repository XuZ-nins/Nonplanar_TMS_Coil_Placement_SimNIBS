%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Xu Zhang @ UConn, Nov. 2022                                                       %
% Determine if each SimNIBS coil is planar (figure-of-8), double-cone, or circular          %
% Dependencies: SimNIBS-4.0                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
allCoils = dir(fullfile(SIMNIBSDIR,'resources\coil_models','Drakaki_BrainStim_2022'));
allCoils_ind = zeros(length(allCoils),1);
for i = 1:length(allCoils)
    if contains(allCoils(i).name,'.ccd')
         allCoils_ind(i) = 1;
    end
end
allCoils(allCoils_ind==0)=[];

for coilnum = 1:length(allCoils)
    tmpcoilname = allCoils(coilnum).name;
    
    tms_opt.fnamecoil = fullfile(SIMNIBSDIR,'resources\coil_models','Drakaki_BrainStim_2022',tmpcoilname);
    
    fileID = fopen(tms_opt.fnamecoil,'r');
    fgetl(fileID);
    dipole_num = str2double(fgetl(fileID));
    coilpos = zeros(dipole_num,3);
    fgetl(fileID);
    for i = 1:dipole_num
        tmpcoilvcs = fgetl(fileID);
        tmpstr = strfind(tmpcoilvcs,' ');
        coilpos(i,1) = str2double(tmpcoilvcs(1:tmpstr(1)-1));
        coilpos(i,2) = str2double(tmpcoilvcs(tmpstr(1)+1:tmpstr(2)-1));
        coilpos(i,3) = str2double(tmpcoilvcs(tmpstr(2)+1:tmpstr(3)-1));
    end
    fclose(fileID);
    coilpos = coilpos*1e3; % Convert to mm
    
    minx0 = min(abs(coilpos_base(:,1)));
    Ix0 = find(abs(abs(coilpos_base(:,1))-minx0)<1);
    maxz0 = max(coilpos_base(Ix0,3));
    maxzl = max(coilpos(coilpos(:,1)>minx0,3));
    maxzr = max(coilpos(coilpos(:,1)<minx0,3));
    rangex = range(coilpos(:,1));
    coilpos_y_nohandle = coilpos(:,2);
    coilpos_y_nohandle(coilpos_y_nohandle>max(abs(coilpos_y_nohandle(coilpos_y_nohandle<0)))) = [];
    coilpos_y_nohandle(coilpos_y_nohandle>max(abs(coilpos_y_nohandle(coilpos_y_nohandle>=0)))) = [];
    rangey = range(coilpos_y_nohandle);
    if rangex>rangey*1.5
        if min(maxzl,maxzr)>maxz0+5
            coilshape = 2; % Double-cone coil
        else
            coilshape = 1; % Figure-of-8 coil
        end
    else
        coilshape = 3; % Circular coil
    end
    switch coilshape
        case 1
            fprintf('%s: Figure-of-8 coil\n',tmpcoilname);
        case 2
            fprintf('%s: Double-cone coil\n',tmpcoilname);
        case 3
            fprintf('%s: Circular coil\n',tmpcoilname);
    end
end
