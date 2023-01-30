function PSFspec = gpuMaskHeight2PSF(t, AP, wavelengths, indices, p_doe, z, nbin, padM, padN)
%GPUMASKHEIGHT2PSF computes spectral PSFs from the height map t masked by AP.
%
%   PSFSPEC = GPUMASKHEIGHT2PSF(T, AP, WAVELENGTHS, INDICES, P_DOE, Z, 
%   NBIN, PADM, PADN) computes the spectral PSFs from the height map t of a  
%   DOE across the WAVELENGTHS range, with the refractive indices of the 
%   INDICES. The DOE is masked by an aperture AP. The DOE resolution is 
%   P_DOE. The propagation distance from DOE substrate to sensor is Z. The 
%   PSFs are actually the Fresnel diffraction patterns for this distance. 
%   The phase profile is zero padded with PADM and PADN before Fresnel 
%   propagration if specified, otherwise no padding is applied. Finally 
%   the PSF intensities are integrated into the sensor pixel area specified 
%   by the ratio between the sensor pixel size and the DOE resolution, 
%   NBIN, which should be a positive integer number. This is the GPU 
%   version of the MaskHeight2PSF functon.

if nargin == 7
    padN = padM;
end

[M, N] = size(t);
K = length(wavelengths);
Mbin = M / nbin;
Nbin = N / nbin;
PSFspec = zeros(Mbin, Nbin, K, 'gpuArray');
for k = 1:K
    phi = 2 * pi * (indices(k) - 1) * t / wavelengths(k);
    Uin = exp(1i*phi) .* AP;
    if nargin > 6
        Uin = padarray(Uin, [padM, padN], 0, 'both');
    end
    Uout = gpuFresnel(Uin, p_doe, p_doe, wavelengths(k), z, 'freq');
    if nargin > 6
        Uout = Uout(padM+1:padM+M, padN+1:padN+N);
    end
    psf = abs(Uout).*2;
    PSFspec(:,:,k) = pixel_integrate(psf, nbin);
end
end