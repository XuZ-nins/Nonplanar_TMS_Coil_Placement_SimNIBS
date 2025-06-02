function [normals_nodes, normals_triangles] = triangle_normals(triangles, nodes, smooth)
% This function calculates the normals of a triangular mesh, both at the 
% nodes (normals_nodes) and at the baricenter of the triangles (normals_triangles).
%
% normals_nodes = triangle_normals(triangles, nodes, smooth);
% [normals_nodes, normals_triangles] = triangle_normals(triangles, nodes, smooth);
%
% Inputs
%   triangles:  A Mx3 list of triangles (faces) of a triangular mesh
%   nodes:      A Nx3 list of vertices of a triangular mesh
%   smooth:     Number of smoothing cycles to perform          
% 
% Outputs
%   normals_nodes: A Nx3 list with the normals of all vertices
%   normals_triangles: A Mx3 list with the normals of all faces
%
% This function is a MATLAB implementation of a Python function within
% mesh_io.py, which is part of the SimNIBS 4.0.1 release
% by Andre Antunes, Guilherme B Saturnino, and Kristoffer H Madsen.
% 
% Written by Xu Zhang at Boston Children's Hospital (June 2024)

    sideA = nodes(double(triangles(:,2)),:)-nodes(triangles(:,1),:);
    sideB = nodes(double(triangles(:,3)),:)-nodes(triangles(:,1),:);
    normals_triangles = cross(sideA, sideB,2);
    normals_triangles = normals_triangles./sqrt(sum(normals_triangles.^2,2));

    % Smooth normals if required
    if smooth > 0
        sideA = nodes(double(triangles(:,2)),:)-nodes(triangles(:,1),:);
        sideB = nodes(double(triangles(:,3)),:)-nodes(triangles(:,1),:);
        normals_triangles = cross(sideA, sideB,2);
    
        normals_nodes = zeros(size(nodes,1),3);
        unique_nodes = unique(triangles);
        for s = 1:smooth+1
            for i = 1:3
                normals_nodes(:,i) = accumarray(reshape(triangles, [],1), ...
                              repmat(normals_triangles(:,i),3,1), ...
                              [size(normals_nodes,1),1]);
            end
            for i = 1:size(triangles,1)
                normals_triangles(i,:) = sum(normals_nodes(triangles(i,:),:),1);
            end
            normals_triangles = normals_triangles./sqrt(sum(normals_triangles.^2,2));
        end
    
        normals_nodes(unique_nodes,:) = normals_nodes(unique_nodes,:)./ ...
            sqrt(sum(normals_nodes(unique_nodes,:).^2,2));
        for i = 1:size(triangles,1)
            normals_triangles(i,:) = sum(normals_nodes(triangles(i,:),:),1);
        end
        normals_triangles = normals_triangles./sqrt(sum(normals_triangles.^2,2));
    end
end
