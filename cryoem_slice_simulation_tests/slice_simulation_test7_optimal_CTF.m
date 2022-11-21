pdbstruct = pdbread( 'helix_40bp.pdb');
pdbstruct = center_pdb( pdbstruct );
% comment out next rotation to see RNA "edge on"
pdbstruct = rotate_pdb( pdbstruct, rotationVectorToMatrix([0,pi/2,0]) );

ice_thickness = 50; % in nm
defocus_vals = [0,0.2,0.5,0.75,1,1.5,2,5,10]*-1; % in um
attenuation = 1; % attenuation by ~10-fold in scattering appears necessary to stay in weak-phase approximation?
all_contrast = [];
for n = 1:length(defocus_vals);
    subplot(5,3,n);
    defocus = defocus_vals(n);
    if n == 1
        [intensity,amplitude] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, ice_thickness );
    else
        [intensity] = simulate_map_from_pdb_FFTbased( amplitude, defocus, ice_thickness );
    end

    all_contrast(n) = max(max(1-intensity));
    title(sprintf('Defocus %5.1f um', defocus))
end
subplot(5,3,10)
[intensity,amplitude] = simulate_map_from_pdb_FFTbased( amplitude, 'phase_plate' );
title( 'Phase plate')
contrast_phase_plate = max(max(1-intensity)); 
subplot(5,3,11)
[intensity,amplitude_fft] = simulate_map_from_pdb_FFTbased( amplitude, 'dark_field' );
title( 'Dark field (note shift in contrast scale)')
contrast_dark_field = max(max(intensity));

subplot(5,1,5)
plot( [0 max(abs(defocus_vals))],contrast_phase_plate*[1 1],'linew',2); hold on
plot( -1 * defocus_vals, all_contrast,'o-','linew',2)
plot( [0 max(abs(defocus_vals))],contrast_dark_field*[1 1],'linew',2); hold off
legend( 'Phase plate','Defocus','Dark field' )
xlabel( '-1 x Defocus (um)');
ylabel( 'Maximum contrast')
title( sprintf('Ice thickness: %5.1f nm',ice_thickness));