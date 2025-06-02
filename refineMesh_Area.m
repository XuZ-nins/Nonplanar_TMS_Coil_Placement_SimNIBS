function [new_vertices, new_faces] = refineMesh_Area(vertices, faces, areaThreshold)
    % Initialize new vertices array with the existing vertices
    new_vertices = vertices;
    edge_midpoint_index = containers.Map('KeyType','char','ValueType','int32');
    
    % Initialize the new faces array
    new_faces = [];

    % Loop over each face
    for i = 1:size(faces, 1)
        face = faces(i, :);
        v1 = vertices(face(1), :);
        v2 = vertices(face(2), :);
        v3 = vertices(face(3), :);
        
        % Calculate the area of the triangle using determinant method
        triArea = 0.5 * abs(v1(1)*(v2(2) - v3(2)) + v2(1)*(v3(2) - v1(2)) + v3(1)*(v1(2) - v2(2)));
        
        if triArea > areaThreshold
            midpoints = zeros(3, 2);

            % Compute midpoints for each edge and check if they are already computed
            for j = 1:3
                vertex_indices = sort(face([j, mod(j, 3) + 1]));
                edge_key = sprintf('%d_%d', vertex_indices);
                
                if edge_midpoint_index.isKey(edge_key)
                    mid_index = edge_midpoint_index(edge_key);
                else
                    midpoint = (vertices(vertex_indices(1), :) + vertices(vertex_indices(2), :)) / 2;
                    new_vertices = [new_vertices; midpoint];
                    mid_index = size(new_vertices, 1);
                    edge_midpoint_index(edge_key) = mid_index;
                end
                midpoints(j, :) = mid_index;
            end
            
            % Create four new triangles from the old triangle
            new_faces = [new_faces; ...
                face(1), midpoints(1), midpoints(3); ...
                midpoints(1), face(2), midpoints(2); ...
                midpoints(3), midpoints(2), face(3); ...
                midpoints(1), midpoints(2), midpoints(3)];
        else
            % If the area does not exceed the threshold, keep the original triangle
            new_faces = [new_faces; face];
        end
    end
end
