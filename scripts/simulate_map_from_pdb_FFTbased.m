function [intensity,amplitude,pixels] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, sigma_blur, max_x, ice_thickness )
% [intensity,amplitude] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, max_x, ice_thickness )
% [intensity,amplitude] = simulate_map_from_pdb_FFTbased( amplitude_fft, max_x, defocus, ice_thickness )
%
%  pdbstruct     = model as read in by pdbread; dimensions assumed to be Å.
%  defocus       = (default -1) defocus, assumed to be in um.
%  sigma_blur    = (default 0 Å) blur by this amount.
%  max_x         = (default 200 Å, 40 nm x 40 nm view) 2D box size will be assumed to be +/- max_x Å.
%  ice_thickness = (default: 20 nm) ice thickness in nm.
%
if ~exist('defocus','var') defocus = -1; end
if ~exist('sigma_blur','var') sigma_blur = 0; end
if ~exist( 'ice_thickness','var') ice_thickness = 20; end;
if ~exist( 'max_x','var') max_x = 200; end;

% pixel size 1 Å
pixel_size = 1e-10;
pixels   = [-max_x:1:max_x]*pixel_size;
N = length(pixels);

[amplitude,lambda]= simulate_amplitude_from_pdb( pdbstruct, max_x, ice_thickness, 2, 0);
amplitude_fft = fft2(amplitude);

if isnumeric(defocus)
    % Let's apply CTF
    midpoint = N/2;
    k = (mod([0:N-1]+midpoint,N)-midpoint)/N;
    [kx,ky] = ndgrid(k,k);

    ctf = exp(-1i*pi*( (defocus*1e-6)*lambda/pixel_size/pixel_size) * (kx.^2+ky.^2));

    amplitude_ctf_fft = amplitude_fft.*ctf;

    % bright field -- add "delta function" at origin.
    N = length(pixels);
    amplitude_ctf_fft(1,1) = amplitude_ctf_fft(1,1) + N*N;
else
    assert( ischar(defocus) )
    if strcmp(defocus,'phase_plate')
        % bright field -- add "delta function", phase shifted by pi/2 at origin.
        amplitude_ctf_fft = amplitude_fft;
        amplitude_ctf_fft(1,1) = amplitude_fft(1,1) - N*N*1i;
    elseif strcmp( defocus, 'dark_field')
        amplitude_ctf_fft = amplitude_fft;
    end
end

amplitude_ctf = ifft2( amplitude_ctf_fft );
intensity = abs(amplitude_ctf).^2;

if (sigma_blur>0.0) intensity = imgaussfilt( intensity, sigma_blur ); end;

%%
% show it!
dark_field = ischar(defocus) & strcmp(defocus,'dark_field');
show_map(intensity, pixels, dark_field);
