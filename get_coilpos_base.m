clear;clc;

addpath(fullfile(SIMNIBSDIR,'matlab_tools'));
addpath(pwd);

cur_dir = pwd;

coil_files = {'MagStim_D70.stl','Deymed_120BFV.stl','MagVenture_Cool-D-B80.stl','MagStim_DCC.stl'};
coil_types = {'D70','120BFV','DB80','DCC'}; % not included in the current study

figure;
for coil_num = 1:length(coil_types)
    coil_type = coil_types{coil_num};
    coil_file = coil_files{coil_num};
    fnamecoil = fullfile(SIMNIBSDIR,'resources','coil_models','Drakaki_BrainStim_2022',coil_file);

    % Read coil geometry
    [F,V] = stlread(fnamecoil);
    	
    mesh = surfaceMesh(V,F);
	fprintf('Initial number of coil vertices: %g\n',size(mesh.Vertices,1));

    [simplified_vertices, simplified_faces] = simplifyMesh_MergeV(mesh.Vertices, mesh.Faces, 1);
    mesh_s1.Vertices = simplified_vertices;
    mesh_s1.Faces = simplified_faces;
    fprintf('Merging close vertices - number of vertices reduced to: %g\n',size(mesh_s1.Vertices,1));
    
    [simplified_vertices, simplified_faces] = simplifyMesh_Area(mesh_s1.Vertices, mesh_s1.Faces, 1);
    mesh_s2.Vertices = simplified_vertices;
    mesh_s2.Faces = simplified_faces;
    fprintf('Merging small faces - number of vertices reduced to: %g\n',size(mesh_s2.Vertices,1));

    [new_vertices, new_faces] = refineMesh_Area(mesh_s2.Vertices, mesh_s2.Faces, 20);
    fprintf('Interpolating large faces - number of vertices increased to: %g\n',size(mesh_s2.Vertices,1));
    
    % coilpos_base = unique(new_vertices,'rows');
    coilpos_base = new_vertices;
    fprintf('Final number of vertices reduced to: %g\n',size(coilpos_base,1));
	coilpos_base = [coilpos_base,ones(size(coilpos_base,1),1)];
    coilpos_face = new_faces;
	% clear F V simplify_flag;

	% Read coil dipoles
	fileID = fopen(replace(fnamecoil,'stl','ccd'),'r');
	fgetl(fileID);
	dipole_num = str2double(fgetl(fileID));
	coilpos_dipole = zeros(dipole_num,3);
	fgetl(fileID);
	for i = 1:dipole_num
	    tmpcoilvcs = fgetl(fileID);
	    tmpstr = strfind(tmpcoilvcs,' ');
	    coilpos_dipole(i,1) = str2double(tmpcoilvcs(1:tmpstr(1)-1));
	    coilpos_dipole(i,2) = str2double(tmpcoilvcs(tmpstr(1)+1:tmpstr(2)-1));
	    coilpos_dipole(i,3) = str2double(tmpcoilvcs(tmpstr(2)+1:tmpstr(3)-1));
	    coilpos_dipole(i,4) = str2double(tmpcoilvcs(tmpstr(3)+1:tmpstr(4)-1));
	    coilpos_dipole(i,5) = str2double(tmpcoilvcs(tmpstr(4)+1:tmpstr(5)-1));
	    coilpos_dipole(i,6) = str2double(tmpcoilvcs(tmpstr(5)+1:end));
    end
    fclose(fileID);
    coilpos_dipole(:,1:3) = coilpos_dipole(:,1:3)*1e3; % Convert to mm

    save(fullfile(cur_dir,strcat('coilpos_base_',coil_type,'.mat')),'coilpos_base','coilpos_face');
    save(fullfile(cur_dir,strcat('coil_dipole_',coil_type,'.mat')),'coilpos_dipole');

    subplot(2,2,coil_num);
    plot3(coilpos_base(:,1),coilpos_base(:,2),coilpos_base(:,3),'.');axis equal;title(coil_type);
end
