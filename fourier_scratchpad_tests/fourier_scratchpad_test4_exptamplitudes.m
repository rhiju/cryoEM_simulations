%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From Su and Coppens, Acta Cryst., (1997). A53, 749-762.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf

% Hydrogen -- fit for 0.0 < s < 2.0 Å^-1
Z = 1; % hydrogen
a = [ 0.43028  0.28537  0.17134 0.09451 0.01725 0.00114];
b = [23.02312 10.20138 51.25444 4.13511 1.35427 0.24269];
assert( abs(Z - sum(a))<0.001 ); 
lambda_mu = 0.37; % screening length
subplot(4,2,1); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'H');

Z = 6; % carbon
a = [ 2.09921 1.80832 1.26159 0.56775 0.26303 0.00010];
b = [13.18997 30.37956 0.69255 0.16831 68.42774 0.44083];
lambda_mu = 0.327;% toggled by hand to get a good fit...
subplot(4,2,2); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'C');

Z = 7; % nitrogen
a = [2.45424 2.15782 1.05782 0.57557 0.44959 0.30480];
b = [18.66694 8.31271 0.46989 42.44646 0.08747 0.47126];
lambda_mu = 0.285;% toggled by hand to get a good fit...
subplot(4,2,3); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'N');

Z = 8; % oxygen
a = [ 2.34752 1.83006 1.61538 1.52402 0.41423 0.26867];
b = [9.69710 18.59876 5.19879 0.32408 39.79099 0.01150];
lambda_mu = 0.255;% toggled by hand to get a good fit...
subplot(4,2,4); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'O');

Z = 15; % phosphorus
a = [6.48197 4.31666  1.73759 1.35793 1.10559 0.00010];
b = [1.89537 27.61455  0.50991  66.28296  0.00010  12.05652];
lambda_mu = 0.305;% toggled by hand to get a good fit... fit is actually not so good.
subplot(4,2,5); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'P');

Z = 11; % sodium
a = [4.16491 2.3807 1.70484 1.59622 0.66291 0.48971];
b = [4.23096 9.48502 0.12559 1.98538 172.13327 82.23091]
lambda_mu = 0.34;% toggled by hand to get a good fit... fit is actually not so good.
subplot(4,2,6); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'Na');

Z = 12; % magnesium
a = [3.90882 2.62159 1.69157 1.52610 1.47907 0.77262];
b = [3.06041 6.12146 0.10357 58.65022 1.56940 125.49980]
lambda_mu = 0.34;% toggled by hand to get a good fit... fit is actually not so good.
subplot(4,2,7); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'Mg');

Z = 17; % chloride
a = [7.13381 6.26972 1.82658 1.62579 0.14431 0.00010];
b = [1.17455 18.57626 0.07869 48.08203 0.07871 23.23894];
lambda_mu = 0.275;% toggled by hand to get a good fit... fit is actually not so good.
assert( abs(Z - sum(a))<0.001 ); 
subplot(4,2,8); plot_electron_scattering_amplitude( Z, a, b, lambda_mu); title( 'Cl');

