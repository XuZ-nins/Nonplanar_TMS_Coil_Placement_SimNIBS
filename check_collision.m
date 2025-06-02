function intersect_count = check_collision(coil_vertices, coil_faces, skin_vertices, skin_faces, ~, dist_thresh)
% Modified collision detection that does not rely on simplifiedIntersect
% Parameters:
%   coil_vertices, coil_faces: Vertices and faces of the coil mesh
%   skin_vertices, skin_faces: Vertices and faces of the skin mesh
%   ~: unused parameter (formerly alphaShape)
%   dist_thresh: Distance threshold for collision detection (default = 4 mm)

% Set default distance threshold if not provided
if nargin < 6 || isempty(dist_thresh)
    dist_thresh = 4; % in mm
end

% Quick rejection based on bounding boxes
coil_bbox = [min(coil_vertices); max(coil_vertices)];
skin_bbox = [min(skin_vertices); max(skin_vertices)];

% Add a margin based on distance threshold
margin = dist_thresh;
if any(coil_bbox(1,:) > skin_bbox(2,:) + margin) || any(coil_bbox(2,:) < skin_bbox(1,:) - margin)
    intersect_count = 0;
    return;
end

% Calculate centroids for additional quick rejection
coil_center = mean(coil_vertices);
skin_center = mean(skin_vertices);

% Calculate rough radius estimates
coil_radius_est = max(vecnorm(coil_vertices - coil_center, 2, 2));
skin_radius_est = max(vecnorm(skin_vertices - skin_center, 2, 2));

% Quick rejection based on centroid distance
center_dist = norm(coil_center - skin_center);
if center_dist > (coil_radius_est + skin_radius_est + dist_thresh)
    intersect_count = 0;
    return;
end

% Simple proximity check between vertices
% Sample a subset of vertices for efficiency
vertex_samples = 50;
coil_vertex_indices = round(linspace(1, size(coil_vertices, 1), min(vertex_samples, size(coil_vertices, 1))));
skin_vertex_indices = round(linspace(1, size(skin_vertices, 1), min(vertex_samples, size(skin_vertices, 1))));

% Check if any vertices are too close to each other
for i = 1:length(coil_vertex_indices)
    coil_point = coil_vertices(coil_vertex_indices(i), :);
    
    for j = 1:length(skin_vertex_indices)
        skin_point = skin_vertices(skin_vertex_indices(j), :);
        if norm(coil_point - skin_point) < dist_thresh * 0.5
            intersect_count = 6;
            return;
        end
    end
end

% No collision found
intersect_count = 0;
end