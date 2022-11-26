pdbstruct = pdbread( 'helix_40bp.pdb');
pdbstruct = center_pdb( pdbstruct );
% comment out next rotation to see RNA "edge on"
pdbstruct = rotate_pdb( pdbstruct, rotationVectorToMatrix([0,pi/2,0]) );

max_x = 200; % in Å
ice_thickness = 20; % in nm
set(gcf,'color','white');
defocus = -2; % in um
sigma_blur = 2.0; % in Angstroms
set(figure(1),'pos',[211   563   508   692]); clf;
subplot(3,1,1)
parfor n = 1:100
    [all_intensity(:,:,n),all_amplitude(:,:,n),pixels] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, sigma_blur, max_x, ice_thickness );
end
title( 'Intensity, single micrograph')

%
subplot(3,1,1);
show_map( all_intensity(:,:,1), pixels );
title( 'Intensity, 1 micrograph')

subplot(3,1,2);
show_map( mean(all_intensity(:,:,1:10),3), pixels );
title( 'Intensity, averaged over 10 micrographs')

subplot(3,1,3);
show_map( mean(all_intensity(:,:,:),3), pixels );
title( sprintf('Intensity, averaged over %d micrographs', size(all_intensity,3)) );


%% Check noise from dosing
set(figure(2),'pos',[211   563   508   692]); clf;
dose = 10; 
all_counts = [];
for n = 1:size(all_intensity,3)
    all_counts(:,:,n) = simulate_counts( all_intensity(:,:,n), dose, pixels);
end

%%
subplot(3,1,1);
show_map( imgaussfilt(all_counts(:,:,1)/dose,5), pixels );
title( 'Intensity, 1 micrograph')

subplot(3,1,2);
show_map( imgaussfilt(mean(all_counts(:,:,1:10),3)/dose,5), pixels );
title( 'Intensity, averaged over 10 micrographs')

subplot(3,1,3);
show_map( imgaussfilt(mean(all_counts(:,:,:),3)/dose,5), pixels );
title( sprintf('Intensity, averaged over %d micrographs', size(all_intensity,3)) );



%% apply CTF correction
set(figure(3),'pos',[611   563   508   692]); clf;
all_intensity_ctf_correct = [];
for n = 1:size(all_intensity,3)
    all_intensity_ctf_correct(:,:,n) = ctf_correct( all_intensity(:,:,n), defocus, pixels);
end


subplot(3,1,1);
clim = [-0.05 0.1];
show_map( all_intensity_ctf_correct(:,:,1), pixels, 0,clim );
title( 'Intensity, 1 micrograph')

subplot(3,1,2);
show_map( mean(all_intensity_ctf_correct(:,:,1:10),3), pixels, 0,clim );
title( 'Intensity, averaged over 10 micrographs')

subplot(3,1,3);
show_map( mean(all_intensity_ctf_correct(:,:,:),3), pixels, 0,clim );
title( sprintf('Intensity, averaged over %d micrographs', size(all_intensity,3)) );






