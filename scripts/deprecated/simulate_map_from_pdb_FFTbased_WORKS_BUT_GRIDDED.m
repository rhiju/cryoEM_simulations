function [amplitude,pixels] = simulate_map_from_pdb( pdbstruct )
% [amplitude,pixels] = simulate_map_from_pdb( pdbstruct )

% 10 nm x 10 nm, 10 nm thickness
% pixel size 1 Å?
pixel_size = 1e-10;
max_x = 60; max_z = 100; % in Angstroms
pixels   = [-max_x:1:max_x]*pixel_size;
pixels_z = [-max_z:1:max_z]*pixel_size;
[X,Y,Z] = ndgrid(pixels,pixels,pixels_z);

%
prefactor = (8.99e9)*(1.60217663e-19); % e/(4 pi epsilon_0) --> V*m

% Prepare grid of scattering potentials from helix.
helix_grid = 0*X;
lambda_mu = 1e-10; % 1 Angstrom -- how well-screened the nucleus charge is.
for n = 1:length(pdbstruct.Model.Atom)
    % later update Z to be scattering factor
    switch pdbstruct.Model.Atom(n).element
        case 'H'
            A = 0.5e-10; % in m
        case 'C'
            A = 2.4e-10;
        case 'N'
            A = 2.2e-10;
        case 'O'
            A = 2.0e-10;
        case 'P'
            A = 5.5e-10;
        otherwise
            printf( 'Unknown! %s\n',pdbstruct.Model.Atom(n).element)
    end
    r = [pdbstruct.Model.Atom(n).X,pdbstruct.Model.Atom(n).Y,pdbstruct.Model.Atom(n).Z];
    helix_grid( round(r(1)+max_x)+1,...
                 round(r(2)+max_x)+1,...
                 round(r(3)+max_z)+1 ...
                 ) = A;
end

E = 300e3; % in eV
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 
amplitude_to_volt_m3 = (4*pi*E*lambda^2);  % converts from scattering amplitude (in m) to Volt-m^3
potential_helix_point = helix_grid * amplitude_to_volt_m3/(pixel_size)^3; % now has units of Volt
potential_helix = potential_helix_point;
%potential_helix = imgaussfilt3(potential_helix_point,1); % Use 1 Å smoothing
max(max(max(potential_helix)))

% sum of potentials
potential = potential_helix;
proj_potential_raw = sum(potential,3)*pixel_size;
fprintf( 'Total projected potential raw (V-m): %f\n',sum(sum(proj_potential_raw)) );
%proj_potential = proj_potential_raw;

% % Now apply filter via Fourier transform
proj_potential_raw_fft = fft2( proj_potential_raw );
N = size(proj_potential_raw_fft,1);
midpoint = N/2;
k = 2*pi*(mod([0:N-1]+midpoint,N)-midpoint)/N;
[kx,ky] = ndgrid(k,k);


% % Gaussian -- QUICK SANITY CHECK THAT I UNDERSTAND 2pi and WAVE VECTOR
%sigma_k = 1;
%formfactor_fft = exp(-sigma_k^2 * (kx.^2 + ky.^2)/2);

% Actual "propagator" for screened Coulomb potential
lambda_mu = 1.0; % 1 Angstrom screening
% Following should be normalized to 1. 
% In the future replace the scattering amplitudes for each atom type
% (H,C,N,O) with the numbers vs. wave vector tabulated based on
% Hartree-Fock!
formfactor_fft = (1/4/pi) * (2*pi)*2./(kx.^2+ky.^2+(1/lambda_mu)^2);
proj_potential_fft = proj_potential_raw_fft .* formfactor_fft;
proj_potential = ifft2( proj_potential_fft );
fprintf( 'Total projected potential via FFT (V-m): %f\n',sum(sum(proj_potential)) );
fprintf( 'Maximum projected potential via  FFT (V-m): %f\n',max(max(proj_potential)) );
fprintf( 'Maximum amplitude via  FFT (rad): %f\n',max(max((pi/lambda/E) * proj_potential_fft)));

% amplitude is related to potential...
amplitude_fft = i * (pi/lambda/E) * proj_potential_fft;

% 

% Let's apply CTF
df = 1e-6; % defocus, in Angstroms?
sigma_k = 5; % size of aperture at lens. Note units are basically radians
propagator_fft = exp(i*pi*(df*lambda/pixel_size/pixel_size) * (kx.^2+ky.^2)) .* exp(- (kx.^2+ky.^2) / 2/ sigma_k.^2);

amplitude_fft = amplitude_fft.*propagator_fft;


% bright field -- add "delta function" at origin.
amplitude_fft(1,1) = amplitude_fft(1,1) + N*N; 

amplitude = ifft2( amplitude_fft );

%
% Phase shift -- SHOULD ACTUALLY APPLY IN FOURIER SPACE!!!
%phase_shift = pi * proj_potential/lambda/E; % need to check factor of pi.
% fprintf( 'Maximum phase shift in 1 Å pixel is: %f rad\n',max(max(phase_shift)));
% fprintf( 'Integrated phase shift in 1 Å pixel is: %f rad * Å^2\n',sum(sum(phase_shift)));

% phase plate or CTF -- SHOULD ACTUALLY APPLY IN FOURIER SPACE!!!
%scalefactor = 10; % Angstrom^2 -- convert above phase shift * Angstrom^2 to just phase shift?
%contrast = scalefactor * (abs( exp(-i * phase_shift/scalefactor ) - 1 + i).^2 - 1);
%contrast = -2 * phase_shift;
%contrast = abs(ifft2( i*(pi/lambda/E) * proj_potential_raw_fft .* propagator_fft )).^2;
intensity = abs(amplitude).^2;


%%
% show it!
subplot(1,1,1);
max_contrast= 2;
imagesc( pixels,pixels,intensity',[0 max_contrast]);
%cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = gray(256);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
title( 'Contrast, ÅxÅ pixels')
axis image


