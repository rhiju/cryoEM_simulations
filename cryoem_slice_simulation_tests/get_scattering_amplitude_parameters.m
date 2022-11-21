function [a,b,Z] = get_scattering_amplitude_parameters( element );
% A = get_scattering_amplitude_parameters( element )
%
% Z. Su and P. Coppens, Acta Cryst.  A53, 749-762 (1997).
%
% Six Gaussian fit
%  (fit for 0.0 < s < 2.0 Å^-1)
% 
% Input
%  element = character with element symbol, e.g., 'C'.
% 
% Output:
% a = scattering amplitude for six Gaussians
% b = coefficient in exponenent for six Gaussians
%     
%

a = zeros(1,6);
b = zeros(1,6);
Z = 0;

switch element
    case 'H'
        Z = 1; % hydrogen
        a = [ 0.43028  0.28537  0.17134 0.09451 0.01725 0.00114];
        b = [23.02312 10.20138 51.25444 4.13511 1.35427 0.24269];
    case 'C'
        Z = 6; % carbon
        a = [ 2.09921 1.80832 1.26159 0.56775 0.26303 0.00010];
        b = [13.18997 30.37956 0.69255 0.16831 68.42774 0.44083];
    case 'N'
        Z = 7; % nitrogen
        a = [2.45424 2.15782 1.05782 0.57557 0.44959 0.30480];
        b = [18.66694 8.31271 0.46989 42.44646 0.08747 0.47126];
    case 'O'
        Z = 8; % oxygen
        a = [ 2.34752 1.83006 1.61538 1.52402 0.41423 0.26867];
        b = [9.69710 18.59876 5.19879 0.32408 39.79099 0.01150];
    case 'P'
        Z = 15; % phosphorus
        a = [6.48197 4.31666  1.73759 1.35793 1.10559 0.00010];
        b = [1.89537 27.61455  0.50991  66.28296  0.00010  12.05652];
    case 'Na'
        Z = 11; % sodium
        a = [4.16491 2.3807 1.70484 1.59622 0.66291 0.48971];
        b = [4.23096 9.48502 0.12559 1.98538 172.13327 82.23091];
    case 'Mg'
        Z = 12; % magnesium
        a = [3.90882 2.62159 1.69157 1.52610 1.47907 0.77262];
        b = [3.06041 6.12146 0.10357 58.65022 1.56940 125.49980];
    case 'Cl'
        Z = 17; % chloride
        a = [7.13381 6.26972 1.82658 1.62579 0.14431 0.00010];
        b = [1.17455 18.57626 0.07869 48.08203 0.07871 23.23894];
    otherwise
        printf( 'Unknown element! %s\n',element);
end
assert( abs(Z - sum(a))<0.001 ); 


