function [simplified_vertices, simplified_faces] = simplifyMesh_MergeV(vertices, faces, distanceThreshold)
    num_vertices = size(vertices, 1);
    merge_map = 1:num_vertices;  % Initially, each vertex maps to itself
    
    % Calculate pairwise distances and merge close vertices
    for i = 1:num_vertices
        for j = i+1:num_vertices
            if norm(vertices(i, :) - vertices(j, :)) < distanceThreshold
                merge_map(j) = merge_map(i);  % Merge vertex j into vertex i
            end
        end
    end
    
    % Update merge_map to its final form
    for i = 1:num_vertices
        merge_map(i) = merge_map(merge_map(i));  % Ensure all mappings are direct
    end

    % Create new vertex array
    new_vertices = zeros(num_vertices, size(vertices, 2));
    count = 0;
    for i = 1:num_vertices
        if merge_map(i) == i  % This vertex is a merged vertex
            count = count + 1;
            new_vertices(count, :) = vertices(i, :);
            merge_map(i) = count;  % Update mapping to new index
        else
            merge_map(i) = merge_map(merge_map(i));  % Update indirect mapping
        end
    end
    new_vertices(count+1:end, :) = [];  % Remove unused space
    
    % Update faces to use new vertex indices
    simplified_faces = faces;
    for i = 1:size(faces, 1)
        for j = 1:3
            simplified_faces(i, j) = merge_map(faces(i, j));
        end
    end
    
    % Remove degenerate triangles (where all vertices are the same)
    simplified_faces = unique(sort(simplified_faces, 2), 'rows');  % Sort to identify duplicates
    degenerate = (simplified_faces(:, 1) == simplified_faces(:, 2)) | ...
                 (simplified_faces(:, 2) == simplified_faces(:, 3));
    simplified_faces(degenerate, :) = [];

    % Return the new vertices and faces
    simplified_vertices = new_vertices;
end
