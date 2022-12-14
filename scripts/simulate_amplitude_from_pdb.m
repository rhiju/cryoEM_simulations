function [amplitude,lambda] = simulate_amplitude_from_pdb( pdbstruct, max_x, ice_thickness, sigma_h2o, SHOW_PLOT )
% [amplitude,lambda] = simulate_amplitude_from_pdb( pdbstruct, max_x, ice_thickness, sigma_h2o, SHOW_PLOT)
%
% Simulate amplitude_fft for 300 keV beam. 
%  Note that total scattering  amplitudes are dependent on atom type as
%   looked up in a Yokemura paper, but the dependence on wavevector is
%   dialed in based on fourier transform of a screened Coulomb potential
%   with screening length fixed at 1 Å.
%
% Inputs
%  pdbstruct     = model as read in by pdbread; dimensions assemed to be Å.
%  max_x         = (default 200 Å) 2D box size will be assumed to be +/- max_x Å.
%  ice_thickness = (default: 20 nm) ice thickness, in nm.
%  sigma_h2o     = (default: 2.0 Å) rms displacement in x/y of h2o beam-movement during dose.
%                   [0.35 Å^2 for each e-/Å^2 of electron dose. So 4 Å^2.
%                   for 10 e-/Å^2. See Mcmullan, Ultramicroscopy, 2015]
%  SHOW_PLOT     = [default 1] show amplitude.
%
% Outputs
%  amplitude = image scattered from model, dimensionless.
%               Note: complex-valued (pure phase shift).
%  lambda = wavelength (in m)
%

if ~exist( 'max_x','var') ice_thickness = 200; end;
if ~exist( 'ice_thickness','var') ice_thickness = 20; end;
if ~exist( 'sigma_h2o','var') sigma_h2o = 2; end;
if ~exist( 'SHOW_PLOT','var') SHOW_PLOT = 1; end;

E = 300e3; % in eV
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 

pixel_size = 1e-10; % assume Angstroms!
pixels   = [-max_x:1:max_x]*pixel_size;

% get ready for Fourier series
N = length(pixels);
midpoint = N/2;
q = 2*pi*(mod([0:N-1]+midpoint,N)-midpoint)/N;
[qx,qy] = ndgrid(q,q);

if ~isstruct( pdbstruct )
    % assume user as inputted amplitude as first argument
    amplitude = pdbstruct;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier synthesis of scattering potentials from model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian -- QUICK SANITY CHECK THAT I UNDERSTAND 2pi and WAVE VECTOR
%sigma = 0.1e-10/pixel_size; 
%formfactor_fft = exp(-sigma^2 * (kx.^2 + ky.^2)/2);
elements = {'H','C','N','O','P','Na','Mg','Cl'};
SCREENED_COULOMB_POTENTIAL = 0; % 0 means use actual scattering amplitudes
scattering_amplitude_for_element = get_scattering_amplitude_for_element( elements, SCREENED_COULOMB_POTENTIAL, qx, qy, pixel_size );

% physical constants needed to convert amplitude (m) to V-m3
h = 6.62607015e-34; m = 9.1093837e-31; e= 1.60217663e-19;
amplitude_to_volt_m3 = h^2/(2*pi*m*e);  % converts from scattering amplitude (in m) to Volt-m^3

max_z = round(ice_thickness * 10/2);
pixels   = [-max_x:1:max_x]*1e-10;
pixels_z = [-max_z:1:max_z]*1e-10;
[X,Y,Z] = ndgrid(pixels,pixels,pixels_z);
model_grid = zeros(length(pixels),length(pixels),length(pixels_z));

% Here we go, sum over atoms in Model(1).
proj_potential_fft = zeros(N,N);
for n = 1:length(pdbstruct.Model(1).Atom)
    x = pdbstruct.Model(1).Atom(n).X+max_x;
    y = pdbstruct.Model(1).Atom(n).Y+max_x;
    z = pdbstruct.Model(1).Atom(n).Z+max_z;
    proj_potential_fft = proj_potential_fft + ...
        amplitude_to_volt_m3/(pixel_size)^2 * exp(-1i*(qx*x+qy*y)) .* ...
        scattering_amplitude_for_element{find(strcmp(elements,pdbstruct.Model(1).Atom(n).element))};
    if round(z)+1>0 && round(z)+1 <= length(pixels_z); model_grid( round(x)+1,round(y)+1,round(z)+1 ) = 1; end;
end

proj_potential_fft = proj_potential_fft ;

% sum of potentials
proj_potential = ifft2(proj_potential_fft);
fprintf( '\nTotal  projected potential (V-nm^3): %f\n',sum(sum(proj_potential))*(pixel_size^2)/(1e-9)^3  );
fprintf( 'Maximum projected potential (V-nm): %f\n',max(max(proj_potential))/1e-9 );
fprintf( 'Maximum amplitude           (rad): %f\n',max(max((pi/lambda/E) * proj_potential)));

% amplitude is related to z-projected potential:
amplitude_model_fft = i * (pi/lambda/E) * proj_potential_fft;


% space_fill grid -- useful for excluding ions and water later.
model_grid_blur = imgaussfilt3( model_grid, 1.4 ); % roll a ball with radius 1.4 Å.
exclude_grid = (model_grid_blur/max(max(max(model_grid_blur))) > 0.6); % 0.4?
fprintf( '\nTotal grid points with model: %d\n', sum(sum(sum(model_grid))));
fprintf( 'Total grid points to exclude water: %d\n',sum(sum(sum(exclude_grid))) );

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form factors used for Na, Cl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coulomb potential screened down to dielectric = 80 at lambda_dielectric.
dielectric_h2o = 80;
coulomb_prefactor = (8.99e9)*(1.60217663e-19); % e/(4 pi epsilon_0) --> V*m
lambda_dielectric = 1e-10/pixel_size; % 1 Angstrom -- how well-screened the nucleus charge is.
scattering_amplitude_coulomb_screened = coulomb_prefactor*(4*pi)./(qx.^2+qy.^2+(1/lambda_dielectric)^2);
% This is the more extended potential that is averaged over water, but
% screened by 1/80.
scattering_amplitude_coulomb_extended = coulomb_prefactor*(4*pi)./(qx.^2+qy.^2+(1/max_x)^2);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phosphate potential and counterion atmosphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_phosphate = 0; n_counterion = 0;

potential_charge = zeros(length(pixels),length(pixels),length(pixels_z));
proj_potential_phosphate_screened_fft = 0*amplitude_model_fft;
proj_potential_phosphate_extended_fft = 0*amplitude_model_fft;
proj_potential_counterion_screened_fft = 0*amplitude_model_fft;
proj_potential_counterion_extended_fft = 0*amplitude_model_fft;
proj_counterion_justatom_fft = 0*amplitude_model_fft;

debye_length = 13e-10/1e-10;

for n = 1:length(pdbstruct.Model(1).Atom)
    % later update Z to be scattering factor
    if strcmp(pdbstruct.Model(1).Atom(n).AtomName ,'OP1') || ...
        strcmp(pdbstruct.Model(1).Atom(n).AtomName ,'OP2')
        Z_charge = -0.5; 
        n_phosphate = n_phosphate+Z_charge;
        x1 = pdbstruct.Model(1).Atom(n).X+max_x;
        y1 = pdbstruct.Model(1).Atom(n).Y+max_x;
        z1 = pdbstruct.Model(1).Atom(n).Z+max_z;
        %R = sqrt((X + max_x*1e-10 - x1*1e-10).^2+(Y + max_x*1e-10 - y1*1e-10).^2+(Z + max_z*1e-10 - z1*1e-10).^2);
        % potential_charge = potential_charge + ...
        %     coulomb_prefactor * Z_charge * (1./R) .*exp(-R/(lambda_dielectric*1e-10));
        proj_potential_phosphate_screened_fft = proj_potential_phosphate_screened_fft + ...
            Z_charge * (1-1/dielectric_h2o) * exp(-1i*(qx*x1+qy*y1)) .* scattering_amplitude_coulomb_screened;
        proj_potential_phosphate_extended_fft = proj_potential_phosphate_extended_fft + ...
            Z_charge/dielectric_h2o * exp(-1i*(qx*x1+qy*y1)) .* scattering_amplitude_coulomb_extended;
    end

    if strcmp(pdbstruct.Model(1).Atom(n).AtomName ,'OP1') 
        Z_charge = 1; 
        dielectric = 1; %80;
        n_counterion = n_counterion+Z_charge;
        % This is a hack Gaussian -- would be much better to use
        % exponential, or even better to use NLPB solution to place
        % counterions.
        x2 = x1 + randn(1)*debye_length; 
        y2 = y1 + randn(1)*debye_length;
        z2 = z1 + randn(1)*debye_length;
        proj_counterion_justatom_fft = proj_counterion_justatom_fft + ...
            amplitude_to_volt_m3/(pixel_size)^2 * scattering_amplitude_for_element{find(strcmp(elements,'Na'))}; 
        proj_potential_counterion_screened_fft = proj_potential_counterion_screened_fft + ...
            Z_charge * (1-1/dielectric_h2o) * exp(-1i*(qx*x2+qy*y2)) .* scattering_amplitude_coulomb_screened;
        proj_potential_counterion_extended_fft = proj_potential_counterion_extended_fft + ...
            Z_charge/dielectric_h2o * exp(-1i*(qx*x2+qy*y2)) .* scattering_amplitude_coulomb_extended;
    end

end
fprintf( '\nTotal charge for phosphate, %f\n',n_phosphate);
fprintf( 'Total charge for counterion, %f\n',n_counterion);

%proj_potential_screened = sum(potential_charge,3)*1e-10;
proj_potential_phosphate_screened = ifft2( proj_potential_phosphate_screened_fft);
fprintf( '\nTotal  projected potential phosphate screened coulomb (V-nm^3): %f\n',sum(sum(proj_potential_phosphate_screened))*(pixel_size^2)/(1e-9)^3 );
fprintf( 'Max    projected potential phosphate screened coulomb (V-nm): %f\n',min(min(proj_potential_phosphate_screened))/1e-9 );
amplitude_potential_phosphate_screened = i * (pi/lambda/E) * proj_potential_phosphate_screened;
amplitude_potential_phosphate_screened_fft = fft2(amplitude_potential_phosphate_screened);

proj_potential_counterion_screened = ifft2( proj_potential_counterion_screened_fft );
fprintf( 'Total  projected potential counterion screened coulomb (V-nm^3): %f\n',sum(sum(proj_potential_counterion_screened))*(pixel_size^2)/(1e-9)^3  );
fprintf( 'Max    projected potential counterion screened coulomb (V-nm): %f\n',max(max(proj_potential_counterion_screened))/1e-9 );
amplitude_potential_counterion_screened = i * (pi/lambda/E) * proj_potential_counterion_screened;
amplitude_potential_counterion_screened_fft = fft2(amplitude_potential_counterion_screened);

proj_potential_phosphate_extended = ifft2( proj_potential_phosphate_extended_fft );
fprintf( '\nTotal  projected potential phosphate extended coulomb (V-nm^3): %f\n',sum(sum(proj_potential_phosphate_extended))*(pixel_size^2)/(1e-9)^3 );
fprintf( 'Max    projected potential phosphate extended coulomb (V-nm): %f\n',min(min(proj_potential_phosphate_extended))/1e-9 );
amplitude_potential_phosphate_extended = i * (pi/lambda/E) * proj_potential_phosphate_extended;

proj_potential_counterion_extended = ifft2( proj_potential_counterion_extended_fft );
fprintf( 'Total  projected potential counterion extended coulomb (V-nm^3): %f\n',sum(sum(proj_potential_counterion_extended))*(pixel_size^2)/(1e-9)^3 );
fprintf( 'Max    projected potential counterion extended coulomb (V-nm): %f\n',min(min(proj_potential_counterion_extended))/1e-9 );

proj_potential_extended = proj_potential_phosphate_extended + proj_potential_counterion_extended;
fprintf( 'Total  projected potential (total) extended coulomb (V-nm^3): %f\n',sum(sum(proj_potential_extended))*(pixel_size^2)/(1e-9)^3  );
fprintf( 'Max    projected potential (total) extended coulomb (V-nm): %f\n',min(min(proj_potential_extended))/1e-9 );

amplitude_potential_counterion_extended = i * (pi/lambda/E) * proj_potential_counterion_extended;

amplitude_potential_screened = amplitude_potential_phosphate_screened + amplitude_potential_counterion_screened;
amplitude_potential_extended = amplitude_potential_phosphate_extended + amplitude_potential_counterion_extended;

amplitude_potential_screened_fft = fft2(amplitude_potential_screened);
amplitude_potential_extended_fft = fft2(amplitude_potential_extended);

proj_counterion_justatom = ifft2( proj_counterion_justatom_fft );
fprintf( '\nTotal  projected potential from Na+ counterions (V-nm^3): %f\n',sum(sum(proj_counterion_justatom))*(pixel_size^2)/(1e-9)^3 );
fprintf( 'Max    projected potential from Na+ counterions (V-nm): %f\n',max(max(proj_counterion_justatom))/1e-9 );
amplitude_counterion_justatom = i * (pi/lambda/E) * proj_counterion_justatom;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Salt -- fill with bulk, add/subtract later
%    based on potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conc_na = 0.1; % in M. So 10 mM
na_per_A3 = conc_na * 6.022e23 * 1e-27;
total_na = round(size(X,1)*size(X,2)*size(X,3)*na_per_A3);
na_grid = 0*X;
for n = 1:total_na
    pos1 = randi(size(X,1));
    pos2 = randi(size(X,2));
    pos3 = randi(size(X,3));
    na_grid(pos1,pos2,pos3) = na_grid(pos1,pos2,pos3) + 1; 
end
fprintf( '\nTotal number of Na+ ions: %d\n',sum(na_grid(:)));
Z_na = 1;
scattering_amplitude_na_total = ...
    amplitude_to_volt_m3/(pixel_size)^2 * scattering_amplitude_for_element{find(strcmp(elements,'Na'))} + ...
    Z_na * (1-1/dielectric_h2o) * scattering_amplitude_coulomb_screened + ...
    Z_na /dielectric_h2o * scattering_amplitude_coulomb_extended;

na_grid_rel = na_grid -  na_per_A3; % This is a constant phase shift.
na_grid_outside_model = na_grid_rel.* (1-exclude_grid);
proj_na_grid = sum( na_grid_outside_model, 3);
proj_potential_na_fft = fft2( proj_na_grid ) .* scattering_amplitude_na_total;

conc_cl = conc_na; % balance Na+ with Cl-
cl_per_A3 = conc_cl * 6.022e23 * 1e-27;
cl_grid = 0*X;
for n = 1:total_na
    pos1 = randi(size(X,1));
    pos2 = randi(size(X,2));
    pos3 = randi(size(X,3));
    cl_grid(pos1,pos2,pos3) = cl_grid(pos1,pos2,pos3) + 1; 
end
fprintf( 'Total number of Cl- ions: %d\n',sum(cl_grid(:)));
Z_cl = -1;
scattering_amplitude_cl_total = ...
    amplitude_to_volt_m3/(pixel_size)^2 * scattering_amplitude_for_element{find(strcmp(elements,'Cl'))} + ...
    Z_cl * (1-1/dielectric_h2o) * scattering_amplitude_coulomb_screened + ...
    Z_cl /dielectric_h2o * scattering_amplitude_coulomb_extended;

cl_grid_rel = cl_grid -  cl_per_A3; % This is a constant phase shift.
cl_grid_outside_model = cl_grid_rel.* (1-exclude_grid);
proj_cl_grid = sum( cl_grid_outside_model, 3);
proj_potential_cl_fft = fft2( proj_cl_grid ) .* scattering_amplitude_cl_total;

proj_potential_salt_fft = proj_potential_na_fft + proj_potential_cl_fft;
proj_potential_salt = ifft2( proj_potential_salt_fft );

amplitude_salt_fft = i * (pi/lambda/E) * proj_potential_salt_fft;
amplitude_salt = ifft2( amplitude_salt_fft );

%proj_potential_extended = proj_potential_phosphate_extended + proj_potential_counterion_extended;
fprintf( '\nTotal  projected potential (bulk NaCl, %5.3f M) (V-nm^3): %f\n',conc_na,sum(sum(proj_potential_salt))*(pixel_size^2)/(1e-9)^3  );
fprintf( 'Max    projected potential (bulk NaCl, %5.3f M) (V-nm): %f\n'  ,conc_na,max(max(proj_potential_salt))/1e-9 );

% TODO 
%  -- add Na only as 1/2 of atmosphere
%  -- decrement Cl (1/2 of atmosphere)
%  -- reduce bulk potential by na+cl contribs.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's go ahead and simulate 50nm thick ice.
% delta functions for potentials of water.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multiply by approximate formfactor for h2o
scattering_amplitude_h2o = ...
    scattering_amplitude_for_element{find(strcmp(elements,'H'))} + ...
    scattering_amplitude_for_element{find(strcmp(elements,'H'))} + ...
    scattering_amplitude_for_element{find(strcmp(elements,'O'))};
% HACK - place hydrogens away from water.
%scattering_amplitude_h2o = ...
%    scattering_amplitude_for_element{find(strcmp(elements,'H'))}.* exp(-i*qx*0.972) + ...
%    scattering_amplitude_for_element{find(strcmp(elements,'H'))}.* exp(-i*qy*0.972) + ...
%    scattering_amplitude_for_element{find(strcmp(elements,'O'))};

% exclude water
h2o_per_A3 = 55 * 6.022e23 * 1e-27;
proj_potential_exclude_water = -1*sum(exclude_grid,3)*h2o_per_A3*amplitude_to_volt_m3/(pixel_size)^2 * scattering_amplitude_h2o(1,1);
proj_potential_exclude_water = 0 * proj_potential_exclude_water;
fprintf( '\nTotal  projected potential exclude water (V-nm^3): %f\n',sum(sum(proj_potential_exclude_water))*(pixel_size^2)/(1e-9)^3 );
fprintf( 'Max    projected potential exclude water (V-nm): %f\n',min(min(proj_potential_exclude_water))/1e-9 );
amplitude_exclude_water = i * (pi/lambda/E) * proj_potential_exclude_water;
amplitude_exclude_water_fft = fft2( amplitude_exclude_water );

h2o_grid = (rand(size(X,1),size(X,2),size(X,3))<h2o_per_A3);
h2o_grid_rel = h2o_grid -  h2o_per_A3; % This is a constant phase shift.
h2o_grid_outside_model = h2o_grid_rel.* (1-exclude_grid);

% summing over z later should give 4 * pi * lambda_mu^2 Z e /(4 pi epsilon0) , with
% units of V-m
potential_h2o = h2o_grid_outside_model * amplitude_to_volt_m3/pixel_size^3 ;
proj_potential_h2o = sum( potential_h2o, 3) * 1e-10;
if (ice_thickness == 0); proj_potential_h2o = 0 * proj_potential_h2o; end;

% smooth by sigma_h2o Angstroms

scattering_amplitude_h2o = scattering_amplitude_h2o .* exp( -sigma_h2o^2 * (qx.^2+qy.^2)/2);
proj_potential_h2o_fft = fft2( proj_potential_h2o ) .* scattering_amplitude_h2o;
    
%     % attempt to simulate all waters. too slow unfortunately!
%     N = length(pixels);
%     n_h2o = round(h2o_per_A3 * N * N* (2*max_z+1))
%     scattering_amplitude_h2o = 2* scattering_amplitude_for_element{find(strcmp(elements,'H'))} + ...
%         scattering_amplitude_for_element{find(strcmp(elements,'O'))};
%     amplitude_h2o_fft = scattering_amplitude_h2o * 0;
%     for n = 1:n_h2o
%         x = rand(2*max_x+1);
%         y = rand(2*max_x+1);
%         amplitude_h2o_fft = amplitude_h2o_fft + ...
%             amplitude_to_volt_m3/(pixel_size)^2 * exp(-i*(qx*x+qy*y)) .* ...
%             scattering_amplitude_h2o;
%     end
%     amplitude_h2o_fft(1,1) = amplitude_h2o_fft(1,1)-...
%         N*N*h2o_per_A3*amplitude_to_volt_m3/(pixel_size)^2 * scattering_amplitude_h2o(1,1);

amplitude_h2o_fft = i * (pi/lambda/E) * proj_potential_h2o_fft;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amplitude_fft = amplitude_model_fft  + amplitude_exclude_water_fft + ...
    amplitude_potential_screened_fft + amplitude_potential_extended_fft + ...
    amplitude_counterion_justatom +  ...
    amplitude_h2o_fft + ...
    amplitude_salt_fft;
amplitude = ifft2( amplitude_fft );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_contrast = 0.1;
if SHOW_PLOT
    subplot(3,3,1); show_image(ifft2(amplitude_model_fft), max_contrast, pixels,'Model')
    subplot(3,3,2); show_image(amplitude_exclude_water, max_contrast, pixels,'Excluded water')
    subplot(3,3,3); show_image(amplitude_potential_phosphate_screened, max_contrast, pixels,'Screened potential from phosphate');
    subplot(3,3,4); show_image(ifft2(amplitude_model_fft+amplitude_exclude_water_fft+amplitude_potential_phosphate_screened_fft), max_contrast, pixels,'Model with excluded water,\newline screened phosphate potential');
    subplot(3,3,5); show_image(ifft2(amplitude_potential_extended_fft), max_contrast, pixels,'Extended potential from phosphate,\newline screened by counterions');
    subplot(3,3,6); show_image(amplitude_counterion_justatom + amplitude_potential_counterion_screened, max_contrast, pixels,'Counterions (with their screened potential)');
    subplot(3,3,7); show_image(ifft2(amplitude_salt_fft), max_contrast, pixels,'Salt (bulk)');
    subplot(3,3,8); show_image(ifft2(amplitude_h2o_fft),max_contrast, pixels,sprintf('%d nm ice outside',ice_thickness))
    subplot(3,3,9); show_image(amplitude, max_contrast, pixels,'Complete image')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function show_image(amplitude, max_contrast, pixels, image_title)
imagesc( pixels/1e-10, pixels/1e-10, -imag(amplitude'), max_contrast*[-1 1]);
axis image
cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
xlabel('x (Å)');ylabel('y (Å)')
title( image_title )
axis( [-60 60 -60 60])
    