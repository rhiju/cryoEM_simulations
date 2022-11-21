function [amplitude_fft,lambda] = simulate_amplitudeFFT_from_pdb( pdbstruct, max_x, ice_thickness )
% [amplitude_fft,kx,ky,lambda] = simulate_amplitudeFFT_from_pdb( pdbstruct, max_x, ice_thickness)
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
%
% Outputs
%  amplitude_fft = 2D FFT of image scattered from model. Points are 2D grid
%                   of wavevectors (2*pi*j/N), with
%                      j = 0, 1, .. N/2, -N/2,... -1
%  lambda = wavelength (in m)
%

if ~exist( 'max_x','var') ice_thickness = 200; end;
if ~exist( 'ice_thickness','var') ice_thickness = 20; end;

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
    % assume user as inputted amplitude_fft as first argument
    amplitude_fft = pdbstruct;
    return;
end

% Fourier synthesis of scattering potentials from helix.
% Gaussian -- QUICK SANITY CHECK THAT I UNDERSTAND 2pi and WAVE VECTOR
%sigma = 0.1e-10/pixel_size; 
%formfactor_fft = exp(-sigma^2 * (kx.^2 + ky.^2)/2);
elements = {'H','C','N','O','P','Na','Mg','Cl'};

SCREENED_COULOMB_POTENTIAL = 0; % 0 means use actual scattering amplitudes
if SCREENED_COULOMB_POTENTIAL
    % Actual "propagator" for screened Coulomb potential
    % Following should be normalized to 1.
    % In the future replace the scattering amplitudes for each atom type
    % (H,C,N,O) with the numbers vs. wave vector tabulated based on
    % Hartree-Fock!
    lambda_mu = 1e-10/pixel_size; % 1 Angstrom -- how well-screened the nucleus charge is.
    formfactor_fft = (1/4/pi/lambda_mu^2) * (4*pi)./(qx.^2+qy.^2+(1/lambda_mu)^2);
    for n = 1:length(elements)
        A = get_scattering_amplitude_at_zero_s( elements{n} );
        scattering_amplitude_for_element{n} = A*formfactor_fft;
    end
else
    % Use Su, Coppens six Gaussian form for scattering amplitude.
    % Go ahead and precompute...
    m = 9.1093837e-31; h = 6.62607015e-34; k_coulomb = 8.9875517923e9; e = 1.60217663e-19;
    scattering_amplitude_prefactor = m*k_coulomb*e*e/(2*h^2) * pixel_size;
    s = sqrt(qx.^2+qy.^2)/(4*pi) + 1e-6;
    for n = 1:length(elements)
        f_optical = 0*s;
        [a,b] = get_scattering_amplitude_parameters( elements{n} );
        for m = 1:length(a)
            f_optical = f_optical + a(m)*exp(-b(m)*s.^2);
        end
        f_electron = scattering_amplitude_prefactor * (sum(a) - f_optical)./s.^2;
        scattering_amplitude_for_element{n} = f_electron * 1e-10; % Su,Coppens output is in Angstroms.
    end
end

% physical constants needed to convert amplitude (m) to V-m3
h = 6.62607015e-34; m = 9.1093837e-31; e= 1.60217663e-19;
amplitude_to_volt_m3 = h^2/(2*pi*m*e);  % converts from scattering amplitude (in m) to Volt-m^3

proj_potential_fft = zeros(N,N);
for n = 1:length(pdbstruct.Model.Atom)
    x = pdbstruct.Model.Atom(n).X+max_x;
    y = pdbstruct.Model.Atom(n).Y+max_x;
    proj_potential_fft = proj_potential_fft + ...
        amplitude_to_volt_m3/(pixel_size)^2 * exp(-i*(qx*x+qy*y)) .* ...
        scattering_amplitude_for_element{find(strcmp(elements,pdbstruct.Model.Atom(n).element))};
end

proj_potential_fft = proj_potential_fft ;

% sum of potentials
proj_potential = ifft2(proj_potential_fft);
fprintf( '\nTotal  projected potential (V-m): %f\n',sum(sum(proj_potential)) );
fprintf( 'Maximum projected potential (V-m): %f\n',max(max(proj_potential)) );
fprintf( 'Maximum amplitude           (rad): %f\n',max(max((pi/lambda/E) * proj_potential_fft)));

% amplitude is related to z-projected potential:
amplitude_model_fft = i * (pi/lambda/E) * proj_potential_fft;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's go ahead and simulate 50nm thick ice.
% delta functions for potentials of water.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_z = round(ice_thickness * 10/2);
ICE_GRID_MODEL = 1;
h2o_per_A3 = 55 * 6.022e23 * 1e-27;
if ICE_GRID_MODEL
    pixels   = [-max_x:1:max_x]*1e-10;
    pixels_z = [-max_z:1:max_z]*1e-10;
    [X,Y,Z] = ndgrid(pixels,pixels,pixels_z);

    h2o_grid = (rand(size(X,1),size(X,2),size(X,3))<h2o_per_A3); % this is a constant phase shift.

    % HEY SHOULD ALSO SIMULATE WATER EXCLUSION FROM MOLECULE!
    %h2o_grid( find(helix_grid>0.5) ) = 0;
    % summing over z later should give 4 * pi * lambda_mu^2 Z e /(4 pi epsilon0) , with
    % units of V-m
    potential_h2o = (h2o_grid - h2o_per_A3) * amplitude_to_volt_m3/(1e-10)^3 ;
    proj_potential_h2o = sum( potential_h2o, 3) * 1e-10;

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

    % smooth by sigma_h2o Angstroms
    sigma_h2o = 2.0;
    scattering_amplitude_h2o = scattering_amplitude_h2o .* exp( -sigma_h2o^2 * (qx.^2+qy.^2)/2);
    proj_potential_h2o_fft = fft2( proj_potential_h2o ) .* scattering_amplitude_h2o;
else
    % attempt to simulate all waters. too slow unfortunately
    N = length(pixels);
    n_h2o = round(h2o_per_A3 * N * N* (2*max_z+1))
    scattering_amplitude_h2o = 2* scattering_amplitude_for_element{find(strcmp(elements,'H'))} + ...
        scattering_amplitude_for_element{find(strcmp(elements,'O'))};
    amplitude_h2o_fft = scattering_amplitude_h2o * 0;
    for n = 1:n_h2o
        x = rand(2*max_x+1);
        y = rand(2*max_x+1);
        amplitude_h2o_fft = amplitude_h2o_fft + ...
            amplitude_to_volt_m3/(pixel_size)^2 * exp(-i*(qx*x+qy*y)) .* ...
            scattering_amplitude_h2o;
    end
    amplitude_h2o_fft(1,1) = amplitude_h2o_fft(1,1)-...
        N*N*h2o_per_A3*amplitude_to_volt_m3/(pixel_size)^2 * scattering_amplitude_h2o(1,1);
end

amplitude_h2o_fft = i * (pi/lambda/E) * proj_potential_h2o_fft;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
amplitude_fft = amplitude_model_fft + amplitude_h2o_fft;
