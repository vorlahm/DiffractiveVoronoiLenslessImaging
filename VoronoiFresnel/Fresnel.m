function Uout = Fresnel(Uin, dx, dy, lambda, z, method)
%FRESNEL propagates wave in z direction with the Fresnel approximation.
%   UOUT = FRESNEL(UIN, DX, DY, LAMBDA, Z, METHOD) simulates the Fresnel
%   diffraction of the filed on the source plane Uin through a distance z
%   to get the output field Uout.
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

Lx = dx * N;
Ly = dy * M;

switch method
    case 'freq'
        fx = -1/(2*dx):1/Lx:1/(2*dx)-1/Lx;
        fy = -1/(2*dy):1/Ly:1/(2*dy)-1/Ly;
        [FX, FY] = meshgrid(fx, fy);
        H = exp(1j*2*pi*z/lambda) .* exp(-1j*pi*lambda*z*(FX.^2+FY.^2));
        H = fftshift(H);
    case 'kernel'
        k = 2*pi/lambda;
        x = dx*(-N/2:N/2-1);
        y = dy*(-M/2:M/2-1);
        [X, Y] = meshgrid(x, y);
        h = exp(1j*2*pi*z/lambda)/(1j*lambda*z) .* exp(1j*k/(2*z)*(X.^2+Y.^2));
        H = fft2(fftshift(h)) * dx*dy;
    otherwise
        error('Unknown simulation method.');
end
Uout = ifftshift(ifft2(H .* fft2(fftshift(Uin))));
end