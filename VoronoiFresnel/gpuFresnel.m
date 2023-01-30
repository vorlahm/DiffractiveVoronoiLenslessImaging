function Uout = gpuFresnel(Uin, dx, dy, lambda, z, method)
%GPUFRESNEL propagates wave in z direction with Fresnel approximation.
%   UOUT = GPUFRESNEL(UIN, DX, DY, LAMBDA, Z, METHOD) simulates the Fresnel
%   diffraction of the filed on the source plane Uin through a distance z
%   to get the output field Uout using gpu acceleration.
%
%       Uin: Optical field on the source plane.
%       dx: Sampling interval in x direction.
%       dy: Sampling interval in y direction.
%       lambda: Wavelength.
%       z: Propagation distance.
%       method: Simulation method used. Two options are available. One is
%           to use the transfer function, and the other is to use with the 
%           impulse response. By default, transfer function is used.
%       Uout: Output field on the observation plane.
%
%   This is the GPU version of the Fresnel function.
%
%   [Reference]
%   Voelz DG. Computational fourier optics: a MATLAB tutorial. 
%             Bellingham, Washington: SPIE press; 2011 Jan.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

if nargin < 6
    method = 'freq';
end

[M, N] = size(Uin);
k = 2*pi/lambda;

Lx = dx * N;
Ly = dy * M;

switch method
    case 'freq'
        fx = gpuArray.colon(-1/(2*dx), 1/Lx, 1/(2*dx)-1/Lx);
        fy = gpuArray.colon(-1/(2*dy), 1/Ly, 1/(2*dy)-1/Ly);
        [FX, FY] = meshgrid(fx, fy);
        H = exp(-1j*pi*lambda*z*(FX.^2+FY.^2));
        H = fftshift(H);
    case 'kernel'
        x = dx*gpuArray.colon(-N/2, 1, N/2-1);
        y = dy*gpuArray.colon(-M/2, 1, M/2-1);
        [X, Y] = meshgrid(x, y);
        h = 1/(1j*lambda*z) * exp(1j*k/(2*z)*(X.^2+Y.^2));
        H = fft2(fftshift(h)) * dx * dy;
    otherwise
        error('Unknown simulation method.');
end
Uout = ifftshift(ifft2(H .* fft2(fftshift(Uin))));
end