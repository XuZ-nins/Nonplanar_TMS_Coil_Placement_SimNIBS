function [min_distance, collision_count] = get_coil_skin_distance(coil_vertices, coil_faces, skin_msh, coil_to_scalp_distance)
% Optimized function that calculates the minimum distance between skin vertices and coil triangles
% and estimates the total number of collisions using simplifiedIntersect
% Returns minimum distance and estimated collision count

% Extract skin vertices and faces
skin_vertices = skin_msh.nodes;
skin_faces = skin_msh.triangles;

% Basic sanity check
if isempty(coil_vertices) || isempty(skin_vertices)
    min_distance = inf;
    collision_count = 0;
    return;
end

% Quick bounding box check
coil_bbox = [min(coil_vertices); max(coil_vertices)];
skin_bbox = [min(skin_vertices); max(skin_vertices)];

% Use conservative margin
margin = coil_to_scalp_distance * 3;
if any(coil_bbox(1,:) > skin_bbox(2,:) + margin) || any(coil_bbox(2,:) < skin_bbox(1,:) - margin)
    % For distant meshes, use bounding box distance as an approximation
    min_distance = max(min(coil_bbox(1,:) - skin_bbox(2,:)), min(skin_bbox(1,:) - coil_bbox(2,:)));
    collision_count = 0;
    return;
end

% Calculate minimum distance between skin vertices and coil triangles
[min_distance, closest_skin_vertices] = calculate_vertex_triangle_distance(skin_vertices, coil_vertices, coil_faces);

% If minimum distance is greater than threshold, no collision
if min_distance > coil_to_scalp_distance
    collision_count = 0;
    return;
end

% Otherwise, detect collisions and estimate count
[collision_count, ~] = estimate_collision_count(coil_vertices, coil_faces, skin_vertices, skin_faces, coil_to_scalp_distance, closest_skin_vertices);
end

function [min_dist, closest_skin_vertices] = calculate_vertex_triangle_distance(skin_vertices, coil_vertices, coil_faces)
% Calculate minimum distance between skin vertices and coil triangles
% Optimized to focus only on potentially close vertices with preallocation
% Returns minimum distance and indices of closest skin vertices

% Initialize with a large value
min_dist = inf;
closest_skin_vertices = [];

% Create KD-tree for the coil vertices
coil_kdtree = createns(coil_vertices);

% Step 1: Find the minimum vertex-to-vertex distance as a rough estimate
[~, initial_dists] = knnsearch(coil_kdtree, skin_vertices, 'K', 1);
[rough_min_dist, min_idx] = min(initial_dists);

% Step 2: Get only the skin vertices that are potentially close to the coil
% Use twice the rough minimum distance as threshold, ensuring we don't miss the true minimum
search_radius = rough_min_dist * 2;
[~, dists] = knnsearch(coil_kdtree, skin_vertices, 'K', 1);
potential_skin_indices = find(dists <= search_radius * 1.5);

% Make sure we include at least some vertices
if isempty(potential_skin_indices)
    % Include at least the closest 10 skin vertices
    [~, sorted_indices] = sort(dists);
    potential_skin_indices = sorted_indices(1:min(10, length(sorted_indices)));
end

% Also always include the vertex with minimum distance from the rough estimate
% Preallocate and use logical indexing to avoid growing arrays
if ~ismember(min_idx, potential_skin_indices)
    potential_skin_indices = [potential_skin_indices; min_idx];
end
potential_skin_indices = unique(potential_skin_indices);

% Step 3: Precompute a mapping from coil vertices to triangles for fast lookup
vertex_to_triangles = cell(size(coil_vertices, 1), 1);

% Estimate maximum number of triangles per vertex for preallocation
% A reasonable estimate is 10-20 triangles per vertex for typical meshes
max_triangles_per_vertex = 20;
for i = 1:size(coil_vertices, 1)
    vertex_to_triangles{i} = zeros(max_triangles_per_vertex, 1);
end

% Count triangles per vertex first
triangle_counts = zeros(size(coil_vertices, 1), 1);
for i = 1:size(coil_faces, 1)
    v1 = coil_faces(i, 1);
    v2 = coil_faces(i, 2);
    v3 = coil_faces(i, 3);
    
    triangle_counts(v1) = triangle_counts(v1) + 1;
    triangle_counts(v2) = triangle_counts(v2) + 1;
    triangle_counts(v3) = triangle_counts(v3) + 1;
end

% Reallocate with exact sizes
for i = 1:size(coil_vertices, 1)
    if triangle_counts(i) > 0
        vertex_to_triangles{i} = zeros(triangle_counts(i), 1);
    else
        vertex_to_triangles{i} = [];
    end
end

% Reset counters for filling arrays
triangle_indices = zeros(size(coil_vertices, 1), 1);
for i = 1:size(coil_faces, 1)
    v1 = coil_faces(i, 1);
    v2 = coil_faces(i, 2);
    v3 = coil_faces(i, 3);
    
    % Add triangle to each vertex's list
    triangle_indices(v1) = triangle_indices(v1) + 1;
    vertex_to_triangles{v1}(triangle_indices(v1)) = i;
    
    triangle_indices(v2) = triangle_indices(v2) + 1;
    vertex_to_triangles{v2}(triangle_indices(v2)) = i;
    
    triangle_indices(v3) = triangle_indices(v3) + 1;
    vertex_to_triangles{v3}(triangle_indices(v3)) = i;
end

% Step 4: For each potentially close skin vertex, find the closest triangle
% Maximum size estimates for preallocation
max_nearest_vertices = 5;
max_nearby_triangles = max_nearest_vertices * max_triangles_per_vertex;
% max_potential_triangles = max_nearby_triangles; 

for i = 1:length(potential_skin_indices)
    skin_idx = potential_skin_indices(i);
    skin_point = skin_vertices(skin_idx, :);
    
    % Find nearest coil vertices
    [nearest_vertices, vertex_dists] = knnsearch(coil_kdtree, skin_point, 'K', max_nearest_vertices);
    
    % Only process if these vertices might be closer than our current minimum
    if min(vertex_dists) > min_dist * 1.5
        continue;
    end
    
    % Collect triangles containing these vertices
    nearby_triangles = zeros(max_nearby_triangles, 1);
    triangle_count = 0;
    
    for j = 1:length(nearest_vertices)
        vertex_idx = nearest_vertices(j);
        vertex_triangles = vertex_to_triangles{vertex_idx};
        
        % Add these triangles to our list
        num_triangles = length(vertex_triangles);
        if num_triangles > 0
            nearby_triangles(triangle_count+1:triangle_count+num_triangles) = vertex_triangles;
            triangle_count = triangle_count + num_triangles;
        end
    end
    
    % Trim to actual size
    if triangle_count > 0
        nearby_triangles = nearby_triangles(1:triangle_count);
        nearby_triangles = unique(nearby_triangles);
    else
        nearby_triangles = [];
    end
    
    % Only check triangles that could potentially be closer than current minimum
    max_triangle_edge = 3.0; % Conservative estimate of maximum triangle edge length
    potential_triangles = zeros(length(nearby_triangles), 1);
    potential_count = 0;
    
    for j = 1:length(nearby_triangles)
        tri_idx = nearby_triangles(j);
        tri = coil_faces(tri_idx, :);
        
        % Calculate triangle centroid
        centroid = (coil_vertices(tri(1), :) + coil_vertices(tri(2), :) + coil_vertices(tri(3), :)) / 3;
        
        % Check if this triangle could possibly be closer than current minimum
        centroid_dist = norm(skin_point - centroid);
        if centroid_dist <= min_dist + max_triangle_edge
            potential_count = potential_count + 1;
            potential_triangles(potential_count) = tri_idx;
        end
    end
    
    % Trim to actual size
    if potential_count > 0
        potential_triangles = potential_triangles(1:potential_count);
    else
        potential_triangles = [];
    end
    
    % Check distance to each potentially close triangle
    for j = 1:length(potential_triangles)
        tri_idx = potential_triangles(j);
        tri = coil_faces(tri_idx, :);
        
        % Get triangle vertices
        v1 = coil_vertices(tri(1), :);
        v2 = coil_vertices(tri(2), :);
        v3 = coil_vertices(tri(3), :);
        
        % Calculate point-triangle distance
        dist = point_triangle_distance(skin_point, v1, v2, v3);
        
        % Update minimum distance if this is closer
        if dist < min_dist
            min_dist = dist;
            closest_skin_vertices = skin_idx;
        end
    end
end
end

function [collision_count, collision_regions] = estimate_collision_count(coil_vertices, coil_faces, skin_vertices, skin_faces, dist_thresh, closest_skin_vertices)
% Estimates the total number of collisions by finding initial collisions 
% and then focusing the search around those regions
% Returns estimated collision count and collision regions

% Step 1: Find initial collision to focus the search
[first_collision, coil_tri_idx, skin_tri_idx] = find_first_collision(coil_vertices, coil_faces, skin_vertices, skin_faces, dist_thresh, closest_skin_vertices);

% If no collision found, return 0
if ~first_collision
    collision_count = 0;
    collision_regions = [];
    return;
end

% Step 2: Define regions around the initial collision site
collision_regions = define_collision_regions(coil_vertices, coil_faces, skin_vertices, skin_faces, coil_tri_idx, skin_tri_idx);

% Step 3: Sample and estimate the total collision count
collision_count = estimate_total_collisions(coil_vertices, coil_faces, skin_vertices, skin_faces, dist_thresh, collision_regions);
end

function [collision_found, coil_tri_idx, skin_tri_idx] = find_first_collision(coil_vertices, coil_faces, skin_vertices, skin_faces, dist_thresh, closest_skin_vertices)
% Finds the first collision between the meshes for efficient estimation
% Focuses search starting from the closest vertex

% Initialize
collision_found = false;
coil_tri_idx = 0;
skin_tri_idx = 0;

% Create KD-tree for coil vertices
coil_kdtree = createns(coil_vertices);

% Start search from the closest skin vertex
if ~isempty(closest_skin_vertices)
    start_indices = closest_skin_vertices;
else
    % If no closest vertex provided, find closest vertices
    [~, dists] = knnsearch(coil_kdtree, skin_vertices, 'K', 1);
    [~, sorted_indices] = sort(dists);
    start_indices = sorted_indices(1:min(10, length(sorted_indices)));
end

% Find triangles around the closest vertices
skin_triangles_to_check = [];
for i = 1:length(start_indices)
    vertex_idx = start_indices(i);
    % Find triangles containing this vertex
    tri_indices = find(skin_faces(:,1) == vertex_idx | ...
                       skin_faces(:,2) == vertex_idx | ...
                       skin_faces(:,3) == vertex_idx);
    skin_triangles_to_check = [skin_triangles_to_check; tri_indices];
end
skin_triangles_to_check = unique(skin_triangles_to_check);

% Find closest coil vertices to these skin vertices
close_coil_indices = [];
for i = 1:length(start_indices)
    [indices, ~] = knnsearch(coil_kdtree, skin_vertices(start_indices(i),:), 'K', 5);
    close_coil_indices = [close_coil_indices; indices'];
end
close_coil_indices = unique(close_coil_indices);

% Find coil triangles to check
coil_triangles_to_check = [];
for i = 1:length(close_coil_indices)
    vertex_idx = close_coil_indices(i);
    % Find triangles containing this vertex
    tri_indices = find(coil_faces(:,1) == vertex_idx | ...
                      coil_faces(:,2) == vertex_idx | ...
                      coil_faces(:,3) == vertex_idx);
    coil_triangles_to_check = [coil_triangles_to_check; tri_indices];
end
coil_triangles_to_check = unique(coil_triangles_to_check);

% Check for intersections
for i = 1:length(coil_triangles_to_check)
    c_face = coil_faces(coil_triangles_to_check(i), :);
    V0 = coil_vertices(c_face(1), :);
    V1 = coil_vertices(c_face(2), :);
    V2 = coil_vertices(c_face(3), :);
    
    for j = 1:length(skin_triangles_to_check)
        s_face = skin_faces(skin_triangles_to_check(j), :);
        U0 = skin_vertices(s_face(1), :);
        U1 = skin_vertices(s_face(2), :);
        U2 = skin_vertices(s_face(3), :);
        
        % Quick distance check
        dists = [
            fast_dist(V0, U0), fast_dist(V0, U1), fast_dist(V0, U2);
            fast_dist(V1, U0), fast_dist(V1, U1), fast_dist(V1, U2);
            fast_dist(V2, U0), fast_dist(V2, U1), fast_dist(V2, U2)
        ];
        
        if min(dists) < dist_thresh * 0.5
            collision_found = true;
            coil_tri_idx = coil_triangles_to_check(i);
            skin_tri_idx = skin_triangles_to_check(j);
            return;
        end
        
        % Check for intersection if distance check passed
        if min(dists) < dist_thresh * 1.5
            try
                if simplifiedIntersect(V0, V1, V2, U0, U1, U2)
                    collision_found = true;
                    coil_tri_idx = coil_triangles_to_check(i);
                    skin_tri_idx = skin_triangles_to_check(j);
                    return;
                end
            catch
                % Skip if intersection test fails
            end
        end
    end
end
end

function regions = define_collision_regions(coil_vertices, coil_faces, skin_vertices, skin_faces, coil_tri_idx, skin_tri_idx)
% Defines regions around the initial collision for focused sampling
% Returns a structure with information about collision regions

% Create region around initial collision
regions = struct();

% Get triangle vertices
c_tri = coil_faces(coil_tri_idx, :);
s_tri = skin_faces(skin_tri_idx, :);

% Calculate centroids
coil_centroid = (coil_vertices(c_tri(1),:) + coil_vertices(c_tri(2),:) + coil_vertices(c_tri(3),:)) / 3;
skin_centroid = (skin_vertices(s_tri(1),:) + skin_vertices(s_tri(2),:) + skin_vertices(s_tri(3),:)) / 3;

% Define region parameters
search_radius = 20; % mm

% Find all vertices within region
coil_vertices_in_region = find_vertices_in_radius(coil_vertices, coil_centroid, search_radius);
skin_vertices_in_region = find_vertices_in_radius(skin_vertices, skin_centroid, search_radius);

% Find triangles in regions
regions.coil_triangles = find_triangles_with_vertices(coil_faces, coil_vertices_in_region);
regions.skin_triangles = find_triangles_with_vertices(skin_faces, skin_vertices_in_region);

% Maximum number of triangles to consider
max_triangles = 500;
if length(regions.coil_triangles) > max_triangles
    regions.coil_triangles = regions.coil_triangles(1:max_triangles);
end
if length(regions.skin_triangles) > max_triangles
    regions.skin_triangles = regions.skin_triangles(1:max_triangles);
end

% Store region centers for potential expansion
regions.coil_center = coil_centroid;
regions.skin_center = skin_centroid;
regions.radius = search_radius;
end

function vertices = find_vertices_in_radius(all_vertices, center, radius)
% Finds all vertices within a given radius of the center point
distances = sqrt(sum((all_vertices - center).^2, 2));
vertices = find(distances <= radius);
end

function count = estimate_total_collisions(coil_vertices, coil_faces, skin_vertices, skin_faces, dist_thresh, regions)
% Estimates total collisions based on sampling in collision regions
% Returns estimated collision count

% Initialize collision counter
collision_counter = 0;

% Preallocate distance array
dists = zeros(9, 1);

% Check all triangle pairs in the region
for i = 1:length(regions.coil_triangles)
    c_face = coil_faces(regions.coil_triangles(i), :);
    V0 = coil_vertices(c_face(1), :);
    V1 = coil_vertices(c_face(2), :);
    V2 = coil_vertices(c_face(3), :);
    
    for j = 1:length(regions.skin_triangles)
        s_face = skin_faces(regions.skin_triangles(j), :);
        U0 = skin_vertices(s_face(1), :);
        U1 = skin_vertices(s_face(2), :);
        U2 = skin_vertices(s_face(3), :);
        
        % Quick distance check
        dists(1) = fast_dist(V0, U0);
        dists(2) = fast_dist(V0, U1);
        dists(3) = fast_dist(V0, U2);
        dists(4) = fast_dist(V1, U0);
        dists(5) = fast_dist(V1, U1);
        dists(6) = fast_dist(V1, U2);
        dists(7) = fast_dist(V2, U0);
        dists(8) = fast_dist(V2, U1);
        dists(9) = fast_dist(V2, U2);
        
        min_dist = min(dists);
        
        % Count as collision if very close
        if min_dist < dist_thresh * 0.5
            collision_counter = collision_counter + 1;
            continue;
        end
        
        % Check for intersection if distance check passed
        if min_dist < dist_thresh * 1.5
            try
                if simplifiedIntersect(V0, V1, V2, U0, U1, U2)
                    collision_counter = collision_counter + 1;
                end
            catch
                % Skip if intersection test fails
            end
        end
    end
    
    % Early stopping if we found a significant number of collisions
    if collision_counter >= 50
        % Scale up to estimate total collisions
        total_scaling = (size(coil_faces, 1) * size(skin_faces, 1)) / (length(regions.coil_triangles) * length(regions.skin_triangles));
        
        % Use conservative estimate - don't scale directly
        count = min(100, floor(collision_counter * sqrt(total_scaling)));
        return;
    end
end

% If only a few collisions found, return exact count
count = collision_counter;
end

function triangles = find_triangles_with_vertices(faces, vertex_indices)
% Find all triangles that contain any of the specified vertices
% Uses a fast vectorized approach with preallocation

% Create lookup array
max_idx = max(faces(:));
if max_idx > 1e7  % For extremely large meshes, use a different approach
    % Estimate maximum number of triangles containing our vertices
    estimated_triangles = min(size(faces, 1), length(vertex_indices) * 10);
    triangles = zeros(estimated_triangles, 1);
    triangle_count = 0;
    
    % Find triangles containing any of the vertices
    for i = 1:size(faces, 1)
        if any(ismember(faces(i,:), vertex_indices))
            triangle_count = triangle_count + 1;
            if triangle_count <= estimated_triangles
                triangles(triangle_count) = i;
            else
                % If we exceed our estimate, extend the array
                triangles = [triangles; i];
            end
        end
    end
    
    % Trim to actual size
    triangles = triangles(1:triangle_count);
else
    % Use vectorized approach for normal sized meshes
    vertex_lookup = false(max_idx, 1);
    vertex_lookup(vertex_indices) = true;
    
    % Find triangles with any matching vertex
    has_v1 = vertex_lookup(faces(:,1));
    has_v2 = vertex_lookup(faces(:,2));
    has_v3 = vertex_lookup(faces(:,3));
    
    triangles = find(has_v1 | has_v2 | has_v3);
end
end

function dist = point_triangle_distance(p, v1, v2, v3)
% Calculate minimum distance from point p to triangle (v1, v2, v3)

% Calculate triangle normal
edge1 = v2 - v1;
edge2 = v3 - v1;
normal = cross(edge1, edge2);
area2 = norm(normal);

% Handle degenerate triangle
if area2 < 1e-10
    % Triangle is degenerate, calculate distance to edges
    dist1 = point_line_distance(p, v1, v2);
    dist2 = point_line_distance(p, v2, v3);
    dist3 = point_line_distance(p, v3, v1);
    dist = min([dist1, dist2, dist3]);
    return;
end

% Normalize the normal
normal = normal / area2;

% Calculate signed distance from point to plane
dist_to_plane = dot(normal, p - v1);

% Project the point onto the plane
p_proj = p - dist_to_plane * normal;

% Use barycentric coordinates to check if projection is inside triangle
v0 = v3 - v1;
v1v = v2 - v1;
v2v = p_proj - v1;

dot00 = dot(v0, v0);
dot01 = dot(v0, v1v);
dot02 = dot(v0, v2v);
dot11 = dot(v1v, v1v);
dot12 = dot(v1v, v2v);

inv_denom = 1 / (dot00 * dot11 - dot01 * dot01);
u = (dot11 * dot02 - dot01 * dot12) * inv_denom;
v = (dot00 * dot12 - dot01 * dot02) * inv_denom;

% Check if point is inside triangle
if (u >= 0) && (v >= 0) && (u + v <= 1)
    % Point projection is inside the triangle
    dist = abs(dist_to_plane);
else
    % Point projection is outside the triangle
    % Calculate distance to edges
    dist1 = point_line_distance(p, v1, v2);
    dist2 = point_line_distance(p, v2, v3);
    dist3 = point_line_distance(p, v3, v1);
    dist = min([dist1, dist2, dist3]);
end
end

function dist = point_line_distance(p, v1, v2)
% Calculate minimum distance from point p to line segment (v1, v2)

% Calculate line direction and length squared
line_dir = v2 - v1;
line_length_sq = dot(line_dir, line_dir);

% Handle degenerate line
if line_length_sq < 1e-10
    dist = norm(p - v1);
    return;
end

% Calculate projection parameter t
t = dot(p - v1, line_dir) / line_length_sq;

% Clamp t to line segment
t = max(0, min(1, t));

% Calculate closest point on line
closest = v1 + t * line_dir;

% Calculate distance to closest point
dist = norm(p - closest);
end

function d = fast_dist(x, y)
% Fast calculation of distance between two 3D points
d = sqrt(sum((x-y).^2));
end