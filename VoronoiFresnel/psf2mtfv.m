function mtfv = psf2mtfv(psf)
%PSF2MTFV calculates the MTF volume (MTFv) metric for input PSF.
%
%   MTFV = PSF2MTFV(PSF) first converts the spectral PSF to panchromatic
%   PSF and then calculates the MTF by Fourier transforming the 
%   panchromatic PSF. MTFV is obtained by summing the MTF values over 
%   all frequency components.
% 
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

PSFpan = sum(psf, 3);
PSFpan = PSFpan / sum(PSFpan(:));
MTFpan = fftshift(abs(fft2(fftshift(PSFpan))));
MTFpan = MTFpan / max(MTFpan(:));
mtfv = squeeze(sum(sum(MTFpan,1),2));
end