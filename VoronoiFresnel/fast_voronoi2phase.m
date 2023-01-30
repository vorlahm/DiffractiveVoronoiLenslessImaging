function lensless_phase = fast_voronoi2phase(centers, vertices, M, N, res, lambda0, z)
%FAST_VORONOI2PHASE generates phase profile from the Voronoi diagram.
%
%   LENSLESS_PHASE = FAST_VORONOI2PHASE(CENTERS, VERTICES, M, N, res, 
%   LAMBDA0, Z) generates the phase profile of the DOE by filling each 
%   Voronoi region with a Fresnel lens phase with its center located at the
%   Voronoi center. 
%   centers: Center coordinates of the Voronoi diagram.
%   vertices: Vertices of the Voronoi diagram.
%   M, N: The size of the design space
%   res: The DOE pixel pitch in um.
%   lambda0: Design wavelength in um.
%   z: Distance from the DOE to the sensor, which is also the focal length 
%      of the corresponding Fresnel lens phase.
%
%   This is a fast version of the voronoi2phase function.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

t    = zeros(M, N);
xmin = res * (-N/2+1);
ymin = res * (-M/2+1);
% generate phase profile from region to region
for j = 1:length(vertices)
    % center coordinates
    xc = centers(1,j);
    yc = centers(2,j);
    appoly = sortvert(vertices{j}); % sort the vertices
    xi = appoly(1,:);
    yi = appoly(2,:);
    ximin = min(xi);
    ximax = max(xi);
    yimin = min(yi);
    yimax = max(yi);
    % local coordinates
    N1 = max(floor((ximin-xmin)/res+1), 1);
    N2 = min(ceil((ximax-xmin)/res+1), N);
    M1 = max(floor((yimin-ymin)/res+1), 1);
    M2 = min(ceil((yimax-ymin)/res+1), M);
    x = xmin + res * ((N1:N2)-1);
    y = ymin + res * ((M1:M2)-1);
    [X, Y] = meshgrid(x, y);
    % sub-aperture
    apmask = poly2mask((xi-x(1))/res, (yi-y(1))/res, M2-M1+1, N2-N1+1);
    temp = -2*pi/lambda0*((X-xc).^2+(Y-yc).^2)/(2*z);
    % sub-phase
    t(M1:M2,N1:N2) = t(M1:M2,N1:N2) + temp .* apmask;
end

lensless_phase = mod(t, 2*pi);
end