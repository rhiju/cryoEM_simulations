function   intensity_ctf_correct = ctf_correct_weiner( all_intensity, defocus_vals, pixels )
%  intensity_ctf_correct = ctf_correct_weiner( intensity, defocus_vals, pixels );
%
% Uses Weiner trick of regularizing by multiplying (not dividing) by CTF in
% Fourier space, then dividing by CTF^2.
%

% PROBABLY SHOULD NOT HARD-CODE THIS!!
lambda = 2.24e-12; % 300 keV electron -->2.24 pm

assert( length(defocus_vals) == size(all_intensity,3));

N = length(pixels);
midpoint = N/2;
k = (mod([0:N-1]+midpoint,N)-midpoint)/N;
[kx,ky] = ndgrid(k,k);

pixel_size = (pixels(2)-pixels(1)); % will be in m, so 1e-10!


intensity_ctf_prod_fft = [];
all_ctf2 = [];

for n = 1:size(all_intensity,3)
    defocus = defocus_vals(n);
    ctf = exp(-1i*pi*( (defocus*1e-6)*lambda/pixel_size/pixel_size) * (kx.^2+ky.^2));
    intensity_fft(:,:,n) = fft2( all_intensity(:,:,n) - 1.0 );
    intensity_ctf_prod_fft(:,:,n) = intensity_fft(:,:,n) .* ctf;
    all_ctf2(:,:,n) = abs(ctf).^2;
end
intensity_ctf_correct_fft = sum(intensity_ctf_prod_fft,3)./sum(all_ctf2,3);

intensity_ctf_correct = real(ifft2( i*intensity_ctf_correct_fft) );



