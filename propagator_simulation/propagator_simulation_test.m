box_size = 4096;
midpoint = box_size/2;
x = [1:box_size]-1;
a = zeros(1,box_size);

% beam?
beam_sigma = 100;
a = exp(-(x-midpoint).^4/2/beam_sigma^4); % + exp(-(x-256).^2/2/10^2); 
%a(midpoint+[-200:200]) = 1;
% object in beam
%a(midpoint+[0:50]) = a(midpoint+[0:50])* 0; % absorber
%a(midpoint+[0:20]) = a(midpoint+[0:20])* i; % phase shifter
source_sigma=5;
x1 = midpoint+40;
a = a.*(1-exp(-(x-x1).^2/2/source_sigma^2));

% double slit
%a(midpoint) = 100;
%a(midpoint+10) = 100;
% source_sigma = 5;
% x1 = midpoint;
% a = a + exp(-(x-x1).^2/2/source_sigma^2);
% x2 = midpoint+40;
% a = a + exp(-(x-x2).^2/2/source_sigma^2);
% x3 = midpoint+50;
% a = a + exp(-(x-x3).^2/2/source_sigma^2);
     
% phase shift  at a point source
%a(midpoint) = a(midpoint)*exp(pi/2*i); 
%a(midpoint+20) = a(midpoint+20) * exp(pi/2*i); 

%a = 0*a + 1;
figure(1);clf;
a_fft = fft(a);
plot( real(a_fft) );
subplot(2,1,1)
b = [];
u = [1:box_size]-1; % not actually in use.
clear i
lambda = 0.05; % wavelength
sigma_prop = 1000; % regularize propagator with a Gaussian.
dz = 500; z_all = 0;
b(:,1) = a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propagate from image plane to lens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:100;
    z = n * dz;
    % need to do 'wraparound' to get propagator contributions at negative x.

    % Fresnel propagator (approximate)
    %       propagator = exp( -i * pi * x.^2 /lambda/z - x.^2 / 2/ sigma_prop.^2) * exp(-2*pi*i*z/lambda) * (i/lambda/z) + ...
    %          exp( -i * pi * (x-box_size).^2 /lambda/z - (x-box_size).^2 / 2/ sigma_prop.^2 ) * exp(-2*pi*i*z/lambda) * (i/lambda/z);

    % Kirchoff integral
    r = sqrt( x.^2 + z^2 );
    propagator1 = i/lambda * (exp( - i * 2 * pi * r/lambda)./r) .* (z./r) .* exp(- x.^2 / 2/ sigma_prop.^2); 
    r = sqrt( (x-box_size).^2 + z^2 ); % wraparound
    propagator2 = i/lambda * (exp( - i * 2 * pi * r/lambda)./r) .* (z./r) .* exp(- (x-box_size).^2 / 2/ sigma_prop.^2); 
    propagator = propagator1 + propagator2;

    propagator_fft = fft(propagator);
    
    % analytical Fresnel
%     theta = (x./z);
%     propagator_fft1 = exp(i*pi*z/lambda * theta.^2) .* exp(- theta.^2 / 2/ sigma_prop.^2);
%     theta = ((x-box_size)./z);
%     propagator_fft2 = exp(i*pi*z/lambda * theta.^2) .* exp(- theta.^2 / 2/ sigma_prop.^2);
%     propagator_fft = propagator_fft1 + propagator_fft2;

    b(:,n+1) = ifft(a_fft.*propagator_fft);
    z_all(n+1) = z;
end
%imagesc(real(b));
imagesc(z_all,x,abs(b).^2,[0 max(abs(b(:,end).^2))]);
%imagesc(z_all,x,real(b),[-0.005 0.005]);
axis image;
hold on; plot([0 max(z_all)],[midpoint midpoint],'r')
colormap(gray(100))

%%%%%%%%%%%%%%%%%
% Apply lens
%%%%%%%%%%%%%%%%%
subplot(2,1,2)
clear i;
cla;
b_lens = b(:,end)';
d_o = z_all(end);
d_f = d_o*(0.8); 
%d_f = 1000*4; 
%d_f = 625;
d_i = 1/(1/d_f-1/d_o);
fprintf( 'Object distance: %d, Focal length: %f, Image distance: %f\n',d_o,d_f,d_i)
sigma_lens = 2000;
b_lens = b_lens .* exp( (i * pi * (x-midpoint).^2/lambda/d_f)) .* exp( -(x-midpoint).^2/2/sigma_lens^2);

%%%%%%%%%%%%%%%%%
% Propagate beyond lens
%%%%%%%%%%%%%%%%%
b_lens_fft = fft(b_lens);
clear c;
for n = 1:2000;
    z = n * dz;
    % need to do 'wraparound' to get propagator contributions at negative x.
%       propagator = exp( -i * pi * x.^2 /lambda/z - x.^2 / 2/ sigma_prop.^2) * exp(-2*pi*i*z/lambda) * (i/lambda/z) + ...
%           exp( -i * pi * (x-box_size).^2 /lambda/z - (x-box_size).^2 / 2/ sigma_prop.^2 ) * exp(-2*pi*i*z/lambda) * (i/lambda/z);

    r = sqrt( x.^2 + z^2 );
    propagator1 = i/lambda * (exp( - i * 2 * pi * r/lambda)./r) .* (z./r) .* exp(- x.^2 / 2/ sigma_prop.^2); 
    r = sqrt( (x-box_size).^2 + z^2 ); % wraparound
    propagator2 = i/lambda * (exp( - i * 2 * pi * r/lambda)./r) .* (z./r) .* exp(- (x-box_size).^2 / 2/ sigma_prop.^2); 
    propagator = propagator1 + propagator2;

    propagator_fft = fft(propagator);

    % analytical Fresnel
%     theta = (x./z);
%     propagator_fft1 = exp(i*pi*z/lambda * theta.^2) .* exp(- theta.^2 / 2/ sigma_prop.^2);
%     theta = ((x-box_size)./z);
%     propagator_fft2 = exp(i*pi*z/lambda * theta.^2) .* exp(- theta.^2 / 2/ sigma_prop.^2);
    %propagator_fft = propagator_fft1 + propagator_fft2;

    c(:,n) = ifft(b_lens_fft.*propagator_fft);
    z_all(n) = z;
end
%imagesc(z_all,x,real(c),[0 0.00005]);
imagesc(z_all,x,abs(c).^2,[0 max(abs(c(:,end).^2))]);
axis image;
hold on; plot([0 max(z_all)],[midpoint midpoint],'r')

%
figure(2)
clf

subplot(5,1,1);
plot(x,abs(a)); 
%hold on; plot(x,real(a));  hold off
title( 'input image' );
xlim([0 box_size])

subplot(5,1,2);
plot(x,abs(b(:,end)));
%hold on; plot(x,real(b(:,end))); hold off
title( 'amplitude right before lens' );
xlim([0 box_size])

subplot(5,1,3);
plot(x,abs(circshift(a_fft,midpoint))); 
%hold on; plot(x,real(circshift(a_fft,midpoint))); hold off
title( 'FFT of input image')
xlim([0 box_size])

subplot(5,1,4);
n = d_f/dz;
plot(x,abs(c(:,n))); 
%hold on; plot(x,real(c(:,n))); hold off
title( 'amplitude at back focal plane (should equal FFT)')
xlim([0 box_size])

subplot(5,1,5);
n = round(d_i/dz);
plot(x,abs(c(:,n))); 
%hold on; plot(x,real(c(:,n))); hold off
title( 'amplitude at image distance')
xlim([0 box_size])
