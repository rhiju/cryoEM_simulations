function   intensity_ctf_correct = ctf_correct( intensity, defocus, pixels);
%  intensity_ctf_correct = ctf_correct( intensity, defocus, pixels);

% PROBABLY SHOULD NOT HARD-CODE THIS!!
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 


N = length(pixels);
midpoint = N/2;
k = (mod([0:N-1]+midpoint,N)-midpoint)/N;
[kx,ky] = ndgrid(k,k);

pixel_size = (pixels(2)-pixels(1)); % will be in m, so 1e-10!
intensity_fft = fft2( intensity );
ctf = exp(-1i*pi*( (defocus*1e-6)*lambda/pixel_size/pixel_size) * (kx.^2+ky.^2));
intensity_ctf_correct = imag(ifft2( intensity_fft ./ ctf ));


