% 1 micron x 1 micron, 100 nm thickness
% pixel size 1 nm?
pixels = [-500:1:500];
[X,Y,Z] = ndgrid(pixels,pixels,[-50:1:50]);

%
prefactor = (8.99e9)*(1.60217663e-19) * 1e9; % e/(4 pi epsilon_0) --> V*nm

% delta functions for potentials of water.
h2o_per_nm3 = 55 * 6.022e23 * 1e-24
h2o_grid = sqrt(h2o_per_nm3)*randn(size(X,1),size(X,2),size(X,3)); %  + h2o_per_nm3  % this is a constant phase shift.
Z_h2o = 10; 
%lambda = 2.24e-3; % FIX! 300 keV electron -->2.24 pm. Convert to nm
lambda_mu = 1e-1; % 1 Angstrom = 0.1 nm

% summing over z later should give 4 * pi * lambda_mu^2 Z e /(4 pi epsilon0) , with
% units of V-nm.
potential = h2o_grid * Z_h2o * (4*pi) * prefactor * (lambda_mu)^2 / 1^3; 

% potential for a dipole
Z_charge = 40; % pretend we have a 20-bp duplex.
dielectric = 1; %80;

x1 =  10.5; y1 = 0; z1 = 0;
x2 = -10.5; y2 = 0; z2 = 0;
potential = potential + ...
    prefactor * (Z_charge/dielectric)  * (  1./sqrt((X-x1).^2+(Y-y1).^2+(Z-z1).^2) - ...
        1./sqrt((X-x2).^2+(Y-y2).^2+(Z-z2).^2) );

%
proj_potential = sum(potential,3);  % in V-nm
max(max(proj_potential))

% show it!
% max_V_nm = 200
% imagesc( proj_potential',max_V_nm*[-1 1]);
% colormap(redwhiteblue(-max_V_nm, max_V_nm))
% colorbar();
% title( 'Projected potential (V * nm)')

%%
% Phase shift
E = 300e3; % in eV
lambda = 2.24e-3; % 300 keV electron -->2.24 pm = 2.24e-3 nm
phase_shift = pi * proj_potential/lambda/E;
% fprintf( 'Maximum phase shift in 1 Å pixel is: %f rad\n',max(max(phase_shift)));
% fprintf( 'Integrated phase shift in 1 Å pixel is: %f rad * Å^2\n',sum(sum(phase_shift)));

contrast = 2 * phase_shift;

% show it!
max_contrast= 2;
imagesc( pixels,pixels,contrast',max_contrast*[-1 1]);
colormap(redwhiteblue(-max_contrast, max_contrast))
colorbar();
title( 'Contrast, nm x nm pixels')
axis image
