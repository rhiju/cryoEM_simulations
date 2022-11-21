function [intensity,amplitude] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, ice_thickness )
% [intensity,amplitude] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, ice_thickness )
% [intensity,amplitude] = simulate_map_from_pdb_FFTbased( amplitude_fft, defocus, ice_thickness )
%
%  pdbstruct     = model as read in by pdbread; dimensions assumed to be Å.
%  defocus       = (default -1) defocus, assumed to be in um.
%  ice_thickness = (default: 20 nm) ice thickness in nm.
%
if ~exist('defocus','var') defocus = -1; end
if ~exist( 'ice_thickness','var') ice_thickness = 20; end;

% pixel size 1 Å
pixel_size = 1e-10;
% 40 nm x 40 nm
max_x = 200;  % in pixels, so 400 Å or 40 nm.
pixels   = [-max_x:1:max_x]*pixel_size;
N = length(pixels);

[amplitude,lambda]= simulate_amplitude_from_pdb( pdbstruct, max_x, ice_thickness);
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

%%
% show it!
cla;
max_contrast= 0.5;
clim = 1+max_contrast*[-1 1];
if ischar(defocus) & strcmp(defocus,'dark_field'); clim = [0 max_contrast^2]; end;
imagesc( pixels/1e-10,pixels/1e-10,intensity',clim);
%cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = gray(256);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
title( 'Contrast, ÅxÅ pixels')
xlabel('x (Å)');ylabel('y (Å)')
axis image


