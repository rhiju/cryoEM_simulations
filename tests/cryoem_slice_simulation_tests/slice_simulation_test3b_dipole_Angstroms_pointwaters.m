% 10 nm x 10 nm, 10 nm thickness
% pixel size 1 Å?
pixels = [-50:1:50]*1e-10;
pixels_z = [-250:1:250]*1e-10;
[X,Y,Z] = ndgrid(pixels,pixels,pixels_z);

%
prefactor = (8.99e9)*(1.60217663e-19); % e/(4 pi epsilon_0) --> V*m

% delta functions for potentials of water.
h2o_per_A3 = 55 * 6.022e23 * 1e-27;
%h2o_grid = sqrt(h2o_per_A3)*randn(size(X,1),size(X,2),size(X,3)); %  + h2o_per_nm3  % this is a constant phase shift.
h2o_grid = (rand(size(X,1),size(X,2),size(X,3))>exp(-h2o_per_A3)) - h2o_per_A3; %  + h2o_per_nm3  % this is a constant phase shift.

%%
pdbstruct = readpdb( 'helix_10bp.pdb');

%%
Z_h2o = 10; 
lambda_mu = 1e-10; % 1 Angstrom

% summing over z later should give 4 * pi * lambda_mu^2 Z e /(4 pi epsilon0) , with
% units of V-m
potential_h2o = h2o_grid * Z_h2o * (4*pi) * prefactor * (lambda_mu)^2/(1e-10)^3; 
potential_h2o = imgaussfilt3(potential_h2o,1); % Use 1 Å smoothing
proj_potential_h2o = sum(potential_h2o,3)*1e-10;  % in V-m
%proj_potential_h2o = imgaussfilt(proj_potential_h2o,1);

% potential for a dipole
Z_charge = 40; % pretend we have a 20-bp duplex.
dielectric = 1; %80;

x1 =  20.5e-10; y1 = 0; z1 = 0;
x2 = -20.5e-10; y2 = 0; z2 = 0;
potential_charge = ...
    prefactor * (Z_charge/dielectric)  * (  1./sqrt((X-x1).^2+(Y-y1).^2+(Z-z1).^2) - ...
    1./sqrt((X-x2).^2+(Y-y2).^2+(Z-z2).^2) );

proj_potential_charge = sum(potential_charge,3)*1e-10;  % in V-m

%
proj_potential = proj_potential_h2o + proj_potential_charge;
max(max(proj_potential))

% show it!
% max_V_nm = 200
% imagesc( proj_potential',max_V_nm*[-1 1]);
% colormap(redwhiteblue(-max_V_nm, max_V_nm))
% colorbar();
% title( 'Projected potential (V * nm)')

%
% Phase shift
E = 300e3; % in eV
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 
phase_shift = pi * proj_potential/lambda/E;
% fprintf( 'Maximum phase shift in 1 Å pixel is: %f rad\n',max(max(phase_shift)));
% fprintf( 'Integrated phase shift in 1 Å pixel is: %f rad * Å^2\n',sum(sum(phase_shift)));

contrast = 2 * phase_shift;

% show it!
figure(1)
max_contrast= 10;
imagesc( pixels,pixels,contrast',max_contrast*[-1 1]);
colormap(redwhiteblue(-max_contrast, max_contrast))
colorbar();
title( 'Contrast, ÅxÅ pixels')
axis image

figure(2)
max_contrast= 2;
contrast_smooth = imgaussfilt(contrast,10);
imagesc( pixels,pixels,contrast_smooth',max_contrast*[-1 1]);
colormap(redwhiteblue(-max_contrast, max_contrast))
colorbar();
title( 'Contrast, smoothed, ÅxÅ pixels')
axis image



