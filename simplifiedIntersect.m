function result = simplifiedIntersect(V0, V1, V2, U0, U1, U2)
% SIMPLIFIEDINTERSECT A simplified version of the triangle intersection test
% This is equivalent to the tri_tri_intersect function from the original code

% Compute plane equation of triangle(V0,V1,V2)
E1 = V1 - V0;
E2 = V2 - V0;
N1 = cross(E1, E2);
d1 = -dot(N1, V0);

% Put U0,U1,U2 into plane equation 1
du0 = dot(N1, U0) + d1;
du1 = dot(N1, U1) + d1;
du2 = dot(N1, U2) + d1;

% Coplanarity check
EPSILON = 1e-6;
if abs(du0) < EPSILON, du0 = 0; end
if abs(du1) < EPSILON, du1 = 0; end
if abs(du2) < EPSILON, du2 = 0; end

du0du1 = du0 * du1;
du0du2 = du0 * du2;

if du0du1 > 0 && du0du2 > 0
    % Same sign on all of them + not equal 0
    result = 0;  % No intersection
    return;
end

% Compute plane of triangle (U0,U1,U2)
E1 = U1 - U0;
E2 = U2 - U0;
N2 = cross(E1, E2);
d2 = -dot(N2, U0);

% Put V0,V1,V2 into plane equation 2
dv0 = dot(N2, V0) + d2;
dv1 = dot(N2, V1) + d2;
dv2 = dot(N2, V2) + d2;

if abs(dv0) < EPSILON, dv0 = 0; end
if abs(dv1) < EPSILON, dv1 = 0; end
if abs(dv2) < EPSILON, dv2 = 0; end

dv0dv1 = dv0 * dv1;
dv0dv2 = dv0 * dv2;

if dv0dv1 > 0 && dv0dv2 > 0
    % Same sign on all of them + not equal 0
    result = 0;  % No intersection
    return;
end

% Compute direction of intersection line
D = cross(N1, N2);

% Compute index to the largest component of D
[~, index] = max(abs(D));

% Simplified projection onto L
vp0 = V0(index);
vp1 = V1(index);
vp2 = V2(index);

up0 = U0(index);
up1 = U1(index);
up2 = U2(index);

% Compute interval for triangle 1
isect1 = computeSimpleIntervals(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2);

% Compute interval for triangle 2
isect2 = computeSimpleIntervals(up0, up1, up2, du0, du1, du2, du0du1, du0du2);

% Check if intervals overlap
isect1 = sort(isect1);
isect2 = sort(isect2);

if isect1(2) < isect2(1) || isect2(2) < isect1(1)
    result = 0;
else
    result = 1;
end
end

function isect = computeSimpleIntervals(vp0, vp1, vp2, d0, d1, d2, d0d1, d0d2)
% Compute intervals for projected triangle
isect = zeros(1,2);
if d0d1 > 0
    % Here we know that d0d2 <= 0
    % That is d0, d1 are on the same side, d2 on the other or on the plane
    isect(1) = vp2 + (vp0 - vp2) * d2 / (d2 - d0);
    isect(2) = vp2 + (vp1 - vp2) * d2 / (d2 - d1);
elseif d0d2 > 0
    % Here we know that d0d1 <= 0
    isect(1) = vp1 + (vp0 - vp1) * d1 / (d1 - d0);
    isect(2) = vp1 + (vp2 - vp1) * d1 / (d1 - d2);
elseif d1*d2 > 0 || d0 ~= 0
    % Here we know that d0d1 <= 0 or that d0 != 0
    isect(1) = vp0 + (vp1 - vp0) * d0 / (d0 - d1);
    isect(2) = vp0 + (vp2 - vp0) * d0 / (d0 - d2);
elseif d1 ~= 0
    isect(1) = vp1 + (vp0 - vp1) * d1 / (d1 - d0);
    isect(2) = vp1 + (vp2 - vp1) * d1 / (d1 - d2);
elseif d2 ~= 0
    isect(1) = vp2 + (vp0 - vp2) * d2 / (d2 - d0);
    isect(2) = vp2 + (vp1 - vp2) * d2 / (d2 - d1);
else
    % Triangles are coplanar
    % This part is handled differently in the main function
    isect = [0, 0];
end
end
