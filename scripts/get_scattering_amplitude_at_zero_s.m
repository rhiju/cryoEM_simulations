function  A = get_scattering_amplitude_at_zero_s( element )
% A = get_scattering_amplitude_at_zero_s( element )
% DEPRECATED: use get_scattering_amplitude_parameters() instead.
% 
% Input
%  element = character with element symbol, e.g., 'C'.
% 
% Output:
% A = scattering amplitude at s=0, in Å. 
%      Read out from Yonekura paper.
%

A = 0;

switch element
    case 'H'
        A = 0.5e-10; % in m
    case 'C'
        A = 2.4e-10;
    case 'N'
        A = 2.2e-10;
    case 'O'
        A = 2.0e-10;
    case 'P'
        A = 5.5e-10;
    otherwise
        fprintf( 'Unknown element! %s\n',element);
end