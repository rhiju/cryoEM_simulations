function plot_electron_scattering_amplitude( Z, a, b, lambda_mu); 

box_size = 1024;
midpoint = box_size/2;
q = 2*pi*(mod([0:box_size-1]+midpoint,box_size)-midpoint)/box_size; % spatial frequency.
s = q/(4*pi);  % q = (2 pi) * 2 sin( theta/2 )/lambda.  s = sin (theta/2) / lambda. 

f_optical = 0*s;
for n = 1:length(a)
    f_optical = f_optical + a(n)*exp(-b(n)*s.^2);
end

pixel_size = 1e-10;
m = 9.1093837e-31; h = 6.62607015e-34; k_coulomb = 8.9875517923e9; e = 1.60217663e-19;
prefactor = m*k_coulomb*e*e/(2*h^2) * pixel_size;
f_electron = prefactor * (sum(a) - f_optical)./s.^2;

f_coulomb_screened = prefactor * Z * (4*pi)^2./(q.^2+(1/lambda_mu)^2);

plot(s,f_electron,'o')
hold on
plot(s,f_coulomb_screened,'o')
hold off
legend( 'Derived from Su, Coppens, 1997',sprintf('Screened coulomb potential (%4.3f Å)',lambda_mu),'Location','SouthEast')
xlabel( 's = sin(\theta_{Bragg})/\lambda');
ylabel( 'Amplitude (Å)')
title( 'Elastic scattering factor')
return;


