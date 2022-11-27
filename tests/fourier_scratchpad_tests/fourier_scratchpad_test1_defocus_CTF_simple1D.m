box_size = 1024;
midpoint = box_size/2;
x = [1:box_size]-1;

%Create image
a = zeros(1,box_size);

% beam?
beam_sigma = 100;
%a = exp(-(x-midpoint).^2/2/beam_sigma^2); % + exp(-(x-256).^2/2/10^2); 
%a(midpoint+[-500:500]) = 1;
a = a+1; % bright field.

% double slit
%a(midpoint) = 100;
%a(midpoint+10) = 100;
source_sigma = 5;
x1 = midpoint;
a = a .*exp( i*exp(-(x-x1).^2/2/source_sigma^2) );
x2 = midpoint+40;
a = a .*exp( 0.5*i*exp(-(x-x2).^2/2/source_sigma^2) );
% x3 = midpoint+50;
% a = a + exp(-(x-x3).^2/2/source_sigma^2);


% Let's apply CTF
df = 100; % defocus
lambda = 0.02; %
z = 10000; % to define angles (should be low).
sigma_k = 1; % there's some sensitivity to this "aperture" parameter
k = x./z/lambda;
max_theta = max(x)/z
propagator_fft1 = exp(i*pi*df*lambda * k.^2) .* exp(- k.^2 / 2/ sigma_k.^2);
k = ((x-box_size)./z/lambda);
propagator_fft2 = exp(i*pi*df*lambda * k.^2) .* exp(- k.^2 / 2/ sigma_k.^2);
propagator_fft = propagator_fft1 + propagator_fft2;

a_fft = fft(a);
a_PSF = ifft( a_fft .* propagator_fft);

clf
plot(x,abs(a)); hold on
plot(x,real(a)); hold on
plot(x,abs(a_PSF).^2 )
xlabel( 'x' )
title( sprintf('CTF-applied bright field image, defocus %f, wavelength %f',df,lambda))
legend( 'Input amplitude (flat!)', 'Real value of input', 'CTF-applied image amplitude')
