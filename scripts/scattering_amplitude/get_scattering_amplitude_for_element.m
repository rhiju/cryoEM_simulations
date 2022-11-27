function scattering_amplitude_for_element = get_scattering_amplitude_for_element(elements, SCREENED_COULOMB_POTENTIAL, qx, qy, pixel_size );
% scattering_amplitude_for_element = get_scattering_amplitude_for_element(elements, SCREENED_COULOMB_POTENTIAL, qx, qy, pixel_size );
%
% OUTPUT:
%  scattering_amplitude_for_element = cell of FFT's of scattering amplitude
%                                       for each element. Units are Å.
% (C) R. Das, 2022

if SCREENED_COULOMB_POTENTIAL
    % Actual "propagator" for screened Coulomb potential
    % Following should be normalized to 1.
    lambda_mu = 1e-10/pixel_size; % 1 Angstrom -- how well-screened the nucleus charge is.
    formfactor_fft = (1/4/pi/lambda_mu^2) * (4*pi)./(qx.^2+qy.^2+(1/lambda_mu)^2);
    for n = 1:length(elements)
        A = get_scattering_amplitude_at_zero_s( elements{n} );
        scattering_amplitude_for_element{n} = A*formfactor_fft;
    end
else
    % Use Su, Coppens six Gaussian form for scattering amplitude, based on
    % relativistic Hartree-Fock.
    % Go ahead and precompute...
    m = 9.1093837e-31; h = 6.62607015e-34; k_coulomb = 8.9875517923e9; e = 1.60217663e-19;
    scattering_amplitude_prefactor = m*k_coulomb*e*e/(2*h^2) * pixel_size;
    s = sqrt(qx.^2+qy.^2)/(4*pi);
    for n = 1:length(elements)
        f_optical = 0*s;
        [a,b] = get_scattering_amplitude_parameters( elements{n} );
        for m = 1:length(a)
            f_optical = f_optical + a(m)*exp(-b(m)*s.^2);
        end
        f_electron = scattering_amplitude_prefactor * (sum(a) - f_optical)./s.^2;
        f_electron(1,1) = scattering_amplitude_prefactor * sum(a.*b); % undefined at s=0, so take limit.
        scattering_amplitude_for_element{n} = f_electron * 1e-10; % Su,Coppens output is in Angstroms.
        % fprintf('%f\n', f_electron(1,1) ) % amplitude at s = 0 !
    end
end
