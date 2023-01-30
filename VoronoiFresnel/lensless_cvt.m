function [MTFv_opt, coords_opt] = lensless_cvt(centers0, lb, ub, maxiter, tol, params, useGPU)
%LENSLESS_CVT finds the optimal coordinates to maximize the MTF volume (MTFv).
%
%   [MTFV_OPT, COORDS_OPT] = LENSLESS_CVT(CENTERS0, LB, UB, MAXITER, TOL,
%   PARAMS, USEGPU) finds the optimal site coordinates to maximize the MTFv
%   metric. 

%   centers0: Initial site coordinates to start the optimization. 
%   lb, ub: Boundaries of the optimization space.
%   maxiter: Maximum number of iterations for the Lloyd algorithm.
%   tol: Tolerance to terminate the Lloyd algorithm. 
%   params: A struc that contains the parameters for the lensless DOE. 
%   usegpu: A flag to indicate if to use GPU.
%
%   The lensless cvt algorithm consists of the following steps:
%       1. Initialize with randomly generated site locations.
%       2. Update the site coordinates with the centroids of each Voronoi
%          resgion (the Lloyd algorithm).
%       3. Calculate the current MTF volume metric (MTFv).
%       4. Update the optimal MTFv and the corresponding coordinates if the
%          current MTFv is higher; otherwise continue.
%       5. Terminate the process once the site location changes are within
%          the tolerance range.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

% parse the lensless DOE parameters
M           = params.M;
N           = params.N;
p_doe       = params.p_doe;
sratio      = params.sratio;
wavelengths = params.wavelengths;
indices     = params.indices;
lambda0     = params.lambda0;
n0          = params.n0;
z           = params.z;
padM        = params.padM;
padN        = params.padN;

% generate CVT coordinates with the Lloyd algorithm
coords = lloyd(centers0, lb, ub, maxiter, tol);

% initialization
MTFv_opt = 0;
coords_opt = [];
% iteration
for i = 1:size(coords, 3)
    % 1. update center coordinates
    centers = coords(:,:,i);
    
    % 2. calculate current MTFv
    [~, ~, vertices] = voronoiPolyhedrons(centers, lb, ub);
    %Lensless_phase = voronoi2phase(centers, vertices, M, N, p_doe, lambda0, z);
    Lensless_phase = fast_voronoi2phase(centers, vertices, M, N, p_doe, lambda0, z);
    
    if useGPU
        t = gpuArray(Lensless_phase * lambda0 / (2 * pi * (n0 - 1)));
        gpuPSFspec = gpuHeight2PSF(gpuArray(t), wavelengths, indices, p_doe, z, sratio, padM, padN);
        PSFspec = gather(gpuPSFspec);
    else
        t = Lensless_phase * lambda0 / (2 * pi * (n0 - 1));
        PSFspec = Height2PSF(t, wavelengths, indices, p_doe, z, sratio, padM, padN);
    end
    
    MTFv = psf2mtfv(PSFspec); % optionally MTFv can be normalized by M*N
    
    % 3. update optimal MTFv and optimal coordinates
    if MTFv > MTFv_opt
        MTFv_opt = MTFv;
        coords_opt = centers;
    end
end
end