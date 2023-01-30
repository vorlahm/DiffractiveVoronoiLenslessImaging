function Iout = pixel_integrate(Iin, ratio)
%PIXEL_INTEGRATE integrates ratio x ratio pixel values in Iin.
%   IOUT = PIXEL_INTEGRATE(IIN, RATIO) integrates the pixel values in Iin
%   in each ratio x ratio subregions. The input image size is [M, N, K], 
%   and output image size is [M/ratio, N/ratio, K]. 
%
%   ratio: A positive integer number and dive the input size with no remainder.
%
%   author: Qiang Fu
%   qiang.fu@kaust.edu.sa
%   2023-01-31

[m, n, k] = size(Iin);
p = m / ratio;
q = n / ratio;
Iout = reshape(Iin, [ratio p ratio q k]);
% Iout = sum(sum(Iout, 1), 3); % naive way, but slow
Iout = reshape(sum(sum(Iout, 1), 3), p, q, k);
end