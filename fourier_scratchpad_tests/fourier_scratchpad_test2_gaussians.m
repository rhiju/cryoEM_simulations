box_size = 1024;
midpoint = box_size/2;
x = [1:box_size]-1;

x0 = 100;
a = zeros(1,box_size);
a(find(x==x0)) = 1;

% Let's apply PSF
lambda = 0.02; %
z = 10000; % to define angles (should be low).
sigma_k = 0.1; % there's some sensitivity to this "aperture" parameter
k = x./z/lambda;
max_theta = max(x)/z
propagator_fft1 =  exp(- k.^2 / 2/ sigma_k.^2);
k = ((x-box_size)./z/lambda);
propagator_fft2 = exp(- k.^2 / 2/ sigma_k.^2);
propagator_fft = propagator_fft1 + propagator_fft2;

a_fft = fft(a);
a_PSF_fft = a_fft .* propagator_fft;
a_PSF = ifft( a_fft .* propagator_fft);
propagator = ifft( propagator_fft );

sigma = 1/sigma_k;
a_PSF_analytical = 1/sqrt(2*pi)/sigma * exp( - (x-x0).^2/sigma^2/2);

u = [0:box_size-1]/box_size;
a_PSF_fft_analytical = exp(i*u*x0) .* propagator_fft;

figure(1)
clf
subplot(2,1,1);
plot( a ); hold on
plot( a_PSF );
plot( propagator )
plot( a_PSF_analytical )
title( 'Real space')
legend('a','a after PSF','propagator','analytical')

subplot(2,1,2);
plot( real(a_fft) ); hold on
plot( real(a_PSF_fft) );
plot( real(propagator_fft) );
plot( real(a_PSF_fft_analytical) )

title( 'Fourier space')
legend('a','a after PSF','propagator','analytical')
ylim([-0.8 1.1])
