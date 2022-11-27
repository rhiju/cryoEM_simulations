%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% let's do coulomb (screened)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box_size = 1024;
midpoint = box_size/2;
x = [0:box_size-1];

figure(2);
x0 = 8; 
lambda_mu = 20; % screening length
[xx,yy] = ndgrid(x,x);
rho = sqrt( (xx-x0).^2 + (yy-midpoint).^2 );
potential = exp(-rho/lambda_mu)./rho;
a_PSF_analytical = sum( potential, 2);

% try to get same answer by applying analytical Fourier transform
%  of potential to a point at x0.
a = zeros(box_size,box_size);
a(find(x==x0),1) = 1;

% don't add, just do wraparound
k = 2*pi*(mod(x+midpoint,box_size)-midpoint)/box_size; % spatial frequency.
[kx,ky] = ndgrid(k,k);
propagator_fft = (4*pi)./(kx.^2+ky.^2+(1/lambda_mu)^2);
%propagator_fft = 0*kx + 1; %test only! delta function propagator!

%a_fft = fft2(a); %this works too!
a_fft = exp(-i*kx*x0);
a_PSF_fft = a_fft .* propagator_fft;
a_PSF = ifft2( a_PSF_fft,'symmetric');
propagator = ifft2( propagator_fft,'symmetric' );


clf
subplot(2,1,1);
plot( x,a(:,1),'linew',2 ); hold on
plot( x,a_PSF(:,1),'linew',2 ); 
plot( x,propagator(:,1),'linew',2 )
plot( x,a_PSF_analytical,'linew',2 )
title( 'Real space')
legend('a','a after PSF','propagator','analytical')
xlim([0 30])

subplot(2,2,3);
imagesc( real(a_fft') )
title( 'input image a, 2D FFT (real component)')

subplot(2,2,4);
imagesc( real(propagator_fft') )
title( 'Propagator 2D FFT (real component')
%plot( real(a_fft) ); hold on
%plot( real(a_PSF_fft) );
%plot( real(propagator_fft),'linew',2 );
%plot( real(a_PSF_fft_analytical) )
%ylim([-0.8 20.1])

