function coords = lloyd(centers, lb, ub, maxiter, tol)
%LLOYD computes the centroidal Voronoi tessellation with initial centers.
%
%   COORDS = LLOYD(CENTERS, LB, UB, MAXITER, TOL) computes the centroidal
%   Voronoi tessellation (CVT) by iteratively updating the generators (or 
%   sites) of each Voronoi region by their centroids. 
%
%   centers: Initial generators. 
%   lb, ub: Boundaries. 
%   max_iter: Maximum number of iterations. 
%   tol: Tolerance to terminate the update.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

% compute the vertices from the centers
[~, ~, vertices] = voronoiPolyhedrons(centers, lb, ub);

coords = [];
coords = cat(3, coords, centers);
for i = 1:maxiter
    % update center locations with centroids
    centroids = zeros(2, length(vertices));
    for j = 1:length(vertices)
        % sort the vertices by angle
        polycoords = sortvert(vertices{j});
        % find the centroids
        [centroidX, centroidY] = centroid(polyshape(polycoords(1,:), polycoords(2,:)));
        % record the coordinates
        centroids(:,j) = [centroidX; centroidY];
    end
    coords = cat(3, coords, centroids);
    
    % update vertices
    [~, ~, vertices] = voronoiPolyhedrons(centroids, lb, ub);
    
    % check convergence
    err = norm(coords(:,:,i+1)-coords(:,:,i), 2)/length(vertices);
    if err < tol
        break;
    end
end
end