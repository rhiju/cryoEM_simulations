pdb_files = {'helix_10bp.pdb','helix_20bp.pdb','helix_40bp.pdb','HIV1_DIS_6bg9.pdb','TRibozyme_7ez0.pdb','tRNAphe_1ehz.pdb','SAM_IV_6ues.pdb','FSE_6xrz.pdb','SL5_FARFAR2_model_c.1.1.pdb'};

max_x = 200; % in Å
ice_thickness = 20; % in nm
set(gcf,'color','white');
defocus = -2.0; % in um
sigma_blur = 2.0; % in Angstroms

defocus = 'phase_plate';

set(figure(1),'pos',[210   418   997   837]); 
clear all_intensity all_amplitude
for n = 1:length(pdb_files)
    subplot(3,3,n);
    pdb_file = ['pdb/',pdb_files{n}]
    pdbstruct = pdbread( pdb_file );
    pdbstruct = center_pdb( pdbstruct );
    % comment out next rotation to see RNA "edge on"
    pdbstruct = rotate_pdb( pdbstruct, rotationVectorToMatrix([0,pi/2,0]) );
    [all_intensity(:,:,n),all_amplitude(:,:,n),pixels] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, sigma_blur, max_x, ice_thickness );
    h =title( pdb_files{n} ); set(h,'interp','none')
end
