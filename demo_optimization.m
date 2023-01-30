%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This demo illustrates the optimization of the lensless phase profile  %
%   of a diffractive optical element (DOE) that maximizes the Modulation  %
%   Transfer Function volume (MTFv).                                      %
%                                                                         %
%   See our paper for the details.                                        %
%                                                                         %
%   Fu Q, Yan DM, Heidrich W. Diffractive lensless imaging with optimized %
%   Voronoi-Fresnel phase. Optics Express. 2022 Dec 5;30(25):45807-23.    %
%                                                                         %
%   author: Qiang Fu                                                      %
%   qiang.fu@kaust.edu.sa                                                 %
%   2023-01-31                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;clc

addpath('./3rdparty');
addpath('./VoronoiFresnel')
% load material data
load('./material/FusedSilica.mat', 'wavelengths', 'indices');

fn = 'demo';

% results may vary if a different seed is used
rng(20);

%% system info %%
lambda0 = 0.55;     % design wavelength,um
n0      = 1.46008;  % refractive index at design wavelength
z       = 2000;     % effective distance from DOE to sensor, um
sratio  = 3;        % sampling ratio = p_sensor / p_doe, should be integer
p_cmos  = 3.45;     % pixel size on sensor
p_doe   = p_cmos / sratio;     % pixel size on DOE, um
% sensor size
% Note: larger sizes would take longer time
pM = 240;
pN = 160;
pK = length(wavelengths);
% DOE size
M = pM * sratio;
N = pN * sratio;
% coordinates
x      = p_doe * (-N/2+1:N/2);
y      = p_doe * (-M/2+1:M/2);
[X, Y] = meshgrid(x, y);
% pad size
padM = 0;
padN = 0;

%% polygon bounds from the Voronoi diagram
lb = [x(1); y(1)];      % lower bound 
ub = [x(end); y(end)];  % upper bound

%% optimization parameters %%
maxiter = 100;
tol     = 0.01*p_doe;
useGPU = true; % slow if running on CPU

%% parameter sweep for number of sites%%
params.M           = M;
params.N           = N;
params.p_doe       = p_doe;
params.sratio      = sratio;
params.wavelengths = wavelengths;
params.indices     = indices;
params.lambda0     = lambda0;
params.n0          = n0;
params.z           = z;
params.padM        = padM;
params.padN        = padN;

% range of the number of Voronoi sites
n_start = 12;
dn      = 1;
n_end   = 50;
n_sites = [n_start:dn:n_end];

MTFVs = zeros(length(n_sites), 1);
for i = 1:length(n_sites)
    fprintf('#%d: number of sites = %d\n', i, n_sites(i));
    centers0 = cat(1, 2 * (rand(1, n_sites(i)) - 0.5) * max(x), ...
        2 * (rand(1, n_sites(i)) - 0.5) * max(y));
    [MTFv_opt, coords_opt] = lensless_cvt(centers0, lb, ub, maxiter, tol, params, useGPU);
    MTFVs(i) = MTFv_opt / (pM * pN);
    save(sprintf('./results/%s_f%dum_%dx%d_%dsites_cvt.mat', fn, z, pM, pN, n_sites(i)), 'MTFv_opt', 'coords_opt');
end

figure('color', 'w');
plot(n_sites, MTFVs, 'linewidth', 2);
xlabel('# of sites', 'fontname', 'Times New Roman', 'fontsize', 14);
ylabel('MTFv', 'fontname', 'Times New Roman', 'fontsize', 14);
title(sprintf('%s for %dx%d region', fn, pM, pN));
set(gca,'FontSize',14,'FontName','Times New Roman');
set(gca,'linewidth',2);