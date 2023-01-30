function lensless_phase = voronoi2phase(centers, vertices, M, N, res, lambda0, z)
%VORONOI2PHASE generates the phase profile from the Voronoi diagram.
%
%   LENSLESS_PHASE = VORONOI2PHASE(CENTERS, VERTICES, M, N, res, LAMBDA0,
%   Z) generates the phase profile of the DOE by filling each Voronoi
%   region with a Fresnel lens phase with its center located at the
%   Voronoi center. 
%   
%   centers: Center coordinates of the Voronoi diagram.
%   vertices: Vertices of the Voronoi diagram.
%   M, N: The size of the design space
%   res: The DOE pixel pitch in um.
%   lambda0: Design wavelength in um.
%   z: Distance from the DOE to the sensor, which is also the focal length 
%      of the corresponding Fresnel lens phase.
%
%   Note: this could be slow when the [M, N] is large. 
%   Use fast_voronoi2phase instead.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

% coordinates
x      = res * (-N/2+1:N/2);
y      = res * (-M/2+1:M/2);
[X, Y] = meshgrid(x, y);
t      = zeros(M, N);
% generate phase profile from region to region
for j = 1:length(vertices)
    % sort the vertices
    appoly = sortvert(vertices{j});
    % convert to positive index
    xi = appoly(1,:)/res + N/2 - 1; 
    yi = appoly(2,:)/res + M/2 - 1;
    % sub-aperture
    apmask = poly2mask(xi, yi, M, N); 
    % sub-phase
    temp = -2*pi/lambda0*((X-centers(1,j)).^2+(Y-centers(2,j)).^2)/(2*z);
    t = t + temp .* apmask;
end
lensless_phase = mod(t, 2*pi);
end