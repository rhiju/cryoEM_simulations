pdbstruct = pdbread( 'pdb/helix_40bp.pdb');
pdbstruct = center_pdb( pdbstruct );
% comment out next rotation to see RNA "edge on"
pdbstruct = rotate_pdb( pdbstruct, rotationVectorToMatrix([0,pi/2,0]) );

max_x = 200; % in Å
ice_thickness = 20; % in nm
set(gcf,'color','white');
defocus = -2; % in um
sigma_blur = 2.0; % in Angstroms
subplot(3,1,1)
parfor n = 1:100
    [all_intensity(:,:,n),all_amplitude(:,:,n),pixels] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus, sigma_blur, max_x, ice_thickness );
end

%
set(figure(1),'pos',[211   563   508   692]); clf;
subplot(3,1,1);
show_map( all_intensity(:,:,1), pixels );
title( 'Intensity, 1 micrograph')

subplot(3,1,2);
show_map( mean(all_intensity(:,:,1:10),3), pixels );
title( 'Intensity, averaged over 10 micrographs')

subplot(3,1,3);
show_map( mean(all_intensity(:,:,:),3), pixels );
title( sprintf('Intensity, averaged over %d micrographs', size(all_intensity,3)) );

%% Plot amplitude directly (avoid CTF)
set(figure(2),'pos',[211   563   508   692]); clf;
clim = [-0.1 0.2];
subplot(3,1,1);
show_map( imag(all_amplitude(:,:,1)), pixels, 0, clim);
title( 'Amplitude, 1 micrograph')

subplot(3,1,2);
show_map( imag(mean(all_amplitude(:,:,1:10),3)), pixels, 0, clim );
title( 'Amplitude, averaged over 10 micrographs')

subplot(3,1,3);
show_map( imag(mean(all_amplitude(:,:,:),3)), pixels, 0, clim );
title( sprintf('Amplitude, averaged over %d micrographs', size(all_intensity,3)) );


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
set(figure(3),'pos',[0   563   508   692]); clf;
all_intensity_ctf_correct = [];
for n = 1:size(all_intensity,3)
    all_intensity_ctf_correct(:,:,n) = ctf_correct( all_intensity(:,:,n), defocus, pixels);
end

subplot(4,1,1);
clim = 0.5*[-0.05 0.1];
show_map( all_intensity_ctf_correct(:,:,1), pixels, 0,clim );
title( 'Intensity CTF-corrected, 1 micrograph')

subplot(4,1,2);
show_map( mean(all_intensity_ctf_correct(:,:,1:10),3), pixels, 0,clim );
title( 'Intensity CTF-corrected, averaged over 10 micrographs')

subplot(4,1,3);
show_map( mean(all_intensity_ctf_correct(:,:,:),3), pixels, 0,clim );
title( sprintf('Intensity CTF-corrected, averaged over %d micrographs', size(all_intensity,3)) );

subplot(4,1,4);
mean_intensity_ctf_correct= ctf_correct( mean(all_intensity,3), defocus,pixels);
show_map(mean_intensity_ctf_correct, pixels, 0, clim);
title( sprintf('Intensity, averaged over %d micrographs then CTF-corrected', size(all_intensity,3)) );

%% What's going on with CTF? Get rid of sigma_blur
max_x_large = 200; % 400;
amplitude_noh2o = simulate_amplitude_from_pdb(pdbstruct, max_x_large, 0 );
%%
[intensity_noh2o,~,pixels_large] = simulate_map_from_pdb_FFTbased( amplitude_noh2o, defocus, 0, max_x_large, 0 );

%%
intensity_noh2o_ctf_correct = ctf_correct( intensity_noh2o, defocus,pixels_large);
clim = 0.5*[-0.05 0.1];
subplot(3,1,1);
show_map(imag(amplitude_noh2o), pixels_large, 0, clim);
title('amplitude, no h2o');
subplot(3,1,2);
show_map(intensity_noh2o, pixels_large, 0, 2*clim+1);
title('intensity contrast, no h2o');
subplot(3,1,3);
show_map( real(intensity_noh2o_ctf_correct), pixels_large, 0, clim);
title('intensity CTF-corrected, no h2o');


