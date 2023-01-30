function PSFspec = Height2PSF(t, wavelengths, indices, p_doe, z, nbin, padM, padN)
%HEIGHT2PSF computes spectral PSFs from the height map t.
%
%   PSFSPEC = HEIGHT2PSF(T, WAVELENGTHS, INDICES, P_DOE, Z, NBIN, PADM, 
%   PADN) computes the spectral PSFs from the height map of a DOE across 
%   the WAVELENGTHS range, with the refractive indices of the substrate 
%   INDICES. 
%   
%   t: Height map of the DOE.
%   wavelengths: Wavelength range to be simulated.
%   indices: Refractive indices of the material that correspond to the
%            wavelength range.
%   p_doe: DOE pixel pitch. 
%   z: Propagation distance from DOE to the sensor. 
%   nbin: The sampling ratio between the DOE pixel pitch and the sensor 
%         pixel pitch. Should be a positive integer number.
%   padM, padN: Number of zeros to pad to the DOE phase profile.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

if nargin == 7
    padN = padM;
end

[M, N] = size(t);
K = length(wavelengths);
Mbin = M / nbin;
Nbin = N / nbin;
PSFspec = zeros(Mbin, Nbin, K);
for k = 1:K
    phi = 2 * pi * (indices(k) - 1) * t / wavelengths(k);
    Uin = exp(1i*phi);
    if nargin > 6
        Uin = padarray(Uin, [padM, padN], 0, 'both');
    end
    Uout = Fresnel(Uin, p_doe, p_doe, wavelengths(k), z, 'freq');
    if nargin > 6
        Uout = Uout(padM+1:padM+M, padN+1:padN+N);
    end
    psf = abs(Uout).*2;
    PSFspec(:,:,k) = pixel_integrate(psf, nbin);
end
end