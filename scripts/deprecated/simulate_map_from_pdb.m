function [phase_shift,pixels] = simulate_map_from_pdb( pdbstruct )
% simulate_map_from_pdb( pdbstruct )

% 10 nm x 10 nm, 10 nm thickness
% pixel size 1 Å?
max_x = 60; max_z = 100; % in Angstroms
pixels   = [-max_x:1:max_x]*1e-10;
pixels_z = [-max_z:1:max_z]*1e-10;
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
            A = 0.5e-10; 
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
potential_helix_point = helix_grid * amplitude_to_volt_m3/(1e-10)^3; 
potential_helix = potential_helix_point;
potential_helix = imgaussfilt3(potential_helix_point,1); % Use 1 Å smoothing

% sum of potentials
potential = potential_helix;
proj_potential = sum(potential,3)*1e-10;
fprintf( 'Total projected potential (V-Angstrom^2): %f\n',sum(sum(proj_potential)) );

% % Now apply filter via Fourier transform
% proj_potential_raw_fft = fft2( proj_potential_raw );
% N = size(proj_potential_raw_fft,1);
% wave_vectors = 2*pi*[1:N]/N; % do I need to normalize by length to get a wave vector? 2*pi?
% [U,V] = meshgrid(wave_vectors,wave_vectors);
% sigma = 10;
% % Gaussian -- QUICK SANITY CHECK THAT I UNDERSTAND 2pi and WAVE VECTOR
% % CONVERSION!
% %
% % OH FORGOT TO WRAP AROUND! ARGH.
% % 
% ctf = exp(-sigma^2 * (U.^2 + V.^2)/2);
% proj_potential = 1000*ifft2( proj_potential_raw_fft .* ctf );


%
% Phase shift -- SHOULD ACTUALLY APPLY IN FOURIER SPACE!!!
phase_shift = pi * proj_potential/lambda/E; % need to check factor of pi.
% fprintf( 'Maximum phase shift in 1 Å pixel is: %f rad\n',max(max(phase_shift)));
% fprintf( 'Integrated phase shift in 1 Å pixel is: %f rad * Å^2\n',sum(sum(phase_shift)));

% phase plate or CTF -- SHOULD ACTUALLY APPLY IN FOURIER SPACE!!!
scalefactor = 10; % Angstrom^2 -- convert above phase shift * Angstrom^2 to just phase shift?
contrast = scalefactor * (abs( exp(-i * phase_shift/scalefactor ) - 1 + i).^2 - 1);

%contrast = -2 * phase_shift;

%%
% show it!
subplot(2,2,1);
max_contrast= 10;
imagesc( pixels,pixels,contrast',max_contrast*[-1 1]);
%cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = gray(256);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
title( 'Contrast, ÅxÅ pixels')
axis image

subplot(2,2,2);
max_contrast= 3;
contrast_smooth = imgaussfilt(contrast,5);
imagesc( pixels,pixels,contrast_smooth',max_contrast*[-1 1]);
colormap(cmap);
colorbar();
title( 'Contrast, smoothed, ÅxÅ pixels')
axis image

