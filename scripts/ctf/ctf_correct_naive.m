function   intensity_ctf_correct = ctf_correct_naive( intensity, defocus, pixels, flat_at_zero );
%  intensity_ctf_correct = ctf_correct( intensity, defocus, pixels, flat_at_zero);
%
%

if ~exist( 'flat_at_zero', 'var' ) flat_at_zero = 0; end;
% PROBABLY SHOULD NOT HARD-CODE THIS!!
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 

N = length(pixels);
midpoint = N/2;
k = (mod([0:N-1]+midpoint,N)-midpoint)/N;
[kx,ky] = ndgrid(k,k);

pixel_size = (pixels(2)-pixels(1)); % will be in m, so 1e-10!
intensity_fft = fft2( intensity - 1.0 );
ctf = exp(-1i*pi*( (defocus*1e-6)*lambda/pixel_size/pixel_size) * (kx.^2+ky.^2));

if flat_at_zero
    k_first_optimum = sqrt((1/2)/(-defocus*1e-6*lambda/pixel_size/pixel_size));
    kk = sqrt(kx.^2+ky.^2);
    idx = find( kk(:) <= k_first_optimum );
    ctf(idx) = exp(1i*pi/2);
end

intensity_ctf_correct = real(ifft2( -i*intensity_fft ./ ctf));
%intensity_ctf_correct = ifft2( -i*intensity_fft ./ ctf,'symmetric');



