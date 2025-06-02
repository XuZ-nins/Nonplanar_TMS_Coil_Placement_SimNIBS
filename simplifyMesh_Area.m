function [simplified_vertices, simplified_faces] = simplifyMesh_Area(vertices, faces, areaThreshold)
    % Calculate areas of all triangles
    areas = zeros(size(faces, 1), 1);
    for i = 1:size(faces, 1)
        v1 = vertices(faces(i, 1), :);
        v2 = vertices(faces(i, 2), :);
        v3 = vertices(faces(i, 3), :);
        areas(i) = 0.5 * abs(v1(1)*(v2(2) - v3(2)) + v2(1)*(v3(2) - v1(2)) + v3(1)*(v1(2) - v2(2)));
    end

    % Find triangles with areas above the threshold
    simplified_faces = faces(areas > areaThreshold, :);

    % Extract the vertices used in the large faces
    unique_vertices_indices = unique(simplified_faces(:));
    simplified_vertices = vertices(unique_vertices_indices, :);

    % Create a map from old vertex indices to new indices
    new_indices_map = zeros(max(unique_vertices_indices), 1);
    new_indices_map(unique_vertices_indices) = 1:length(unique_vertices_indices);

    % Update face indices to correspond to the new vertices array
    for i = 1:size(simplified_faces, 1)
        simplified_faces(i, :) = new_indices_map(simplified_faces(i, :));
    end
end

