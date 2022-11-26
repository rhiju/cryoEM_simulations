pdbstruct = pdbread( 'pdb/helix_40bp.pdb');
pdbstruct = pdbread('pdb/TRibozyme_7ez0.pdb');
pdbstruct = center_pdb( pdbstruct );
% comment out next rotation to see RNA "edge on"
pdbstruct = rotate_pdb( pdbstruct, rotationVectorToMatrix([0,pi/2,0]) );

max_x = 200; % in Å
ice_thickness = 20; % in nm
figure(2); amplitude = simulate_amplitude_from_pdb(pdbstruct, max_x, ice_thickness );
set(gcf,'color','white');

figure(1)
set(gcf,'color','white');
sigma_blur = 2.0; % in Angstroms
all_contrast = [];
defocus_vals = [0,0.2,0.5,0.75,1,1.5,2,5,10]*-1; % in um
for n = 1:length(defocus_vals);
    subplot(5,3,n);
    defocus = defocus_vals(n);
    intensity = simulate_map_from_pdb_FFTbased( amplitude, defocus, sigma_blur, max_x, ice_thickness );
    all_contrast(n) = max(max(1-intensity));
    title(sprintf('Defocus %5.1f um', defocus))
end
subplot(5,3,10)
[intensity,amplitude] = simulate_map_from_pdb_FFTbased( amplitude, 'phase_plate', sigma_blur );
title( 'Phase plate')
contrast_phase_plate = max(max(1-intensity)); 
subplot(5,3,11)
[intensity,amplitude_fft] = simulate_map_from_pdb_FFTbased( amplitude, 'dark_field', sigma_blur );
title( 'Dark field\newline(note shift in contrast scale)')
contrast_dark_field = max(max(intensity));

subplot(5,1,5)
plot( [0 max(abs(defocus_vals))],contrast_phase_plate*[1 1],'linew',2); hold on
plot( -1 * defocus_vals, all_contrast,'o-','linew',2)
plot( [0 max(abs(defocus_vals))],contrast_dark_field*[1 1],'linew',2); hold off
legend( 'Phase plate','Defocus','Dark field' )
xlabel( '-1 x Defocus (um)');
ylabel( 'Maximum contrast')
title( sprintf('Ice thickness: %5.1f nm',ice_thickness));
