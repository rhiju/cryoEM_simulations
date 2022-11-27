pdbstruct = pdbread( 'pdb/helix_40bp.pdb');
pdbstruct = center_pdb( pdbstruct );
% comment out next rotation to see RNA "edge on"
pdbstruct = rotate_pdb( pdbstruct, rotationVectorToMatrix([0,pi/2,0]) );

max_x = 200; % in Å
ice_thickness = 20; % in nm
set(gcf,'color','white');
%defocus = -2; % in um
sigma_blur = 0.0; % in Angstroms
subplot(3,1,1)
defocus_vals = -1-[0.01:0.01:1];
parfor n = 1:100
    [all_intensity(:,:,n),all_amplitude(:,:,n),pixels] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus_vals(n), sigma_blur, max_x, ice_thickness );
end

% Need to assign pixels, which we can't get out of parfor.
n = 1; [all_intensity(:,:,n),all_amplitude(:,:,n),pixels] = simulate_map_from_pdb_FFTbased( pdbstruct, defocus_vals(n), sigma_blur, max_x, ice_thickness );

%
set(figure(1),'pos',[211   563   508   692]); clf;
subplot(3,1,1);
show_map( all_intensity(:,:,1), pixels );
title( sprintf('Intensity, 1 micrograph, defocus %5.2f um',defocus_vals(1)))

subplot(3,1,2);
show_map( mean(all_intensity(:,:,1:10),3), pixels );
title( sprintf('Intensity, averaged over 10 micrographs, defocus %5.2f to %5.2f um', defocus_vals(1),defocus_vals(10)) )

subplot(3,1,3);
show_map( mean(all_intensity(:,:,:),3), pixels );
title( sprintf('Intensity, averaged over %d micrographs, defocus %5.2f to %5.2f um', size(all_intensity,3),defocus_vals(1),defocus_vals(end)) )

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
set(figure(3),'pos',[211   563   508   692]); clf;
dose = 10; 
all_counts = [];
for n = 1:size(all_intensity,3)
    all_counts(:,:,n) = simulate_counts( all_intensity(:,:,n), dose, pixels);
end

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
set(figure(4),'pos',[0   563   508   692]); clf;
all_intensity_ctf_correct = [];
for n = 1:size(all_intensity,3)
    all_intensity_ctf_correct(:,:,n) = ctf_correct_naive( all_intensity(:,:,n), defocus_vals(n), pixels);
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
mean_intensity_ctf_correct= ctf_correct_naive( mean(all_intensity,3), mean(defocus_vals),pixels);
show_map(mean_intensity_ctf_correct, pixels, 0, clim);
title( sprintf('Intensity, averaged over %d micrographs then CTF-corrected', size(all_intensity,3)) );


%% Better CTF corrections
%apply CTF correction "flatten" to -1 at k less than first optimum. 
set(figure(6),'pos',[0   563   508   692]); clf;

mean_intensity_ctf_correct_naive =  mean(all_intensity_ctf_correct(:,:,:),3);

all_intensity_ctf_correct_flatten = [];
for n = 1:size(all_intensity,3)
    all_intensity_ctf_correct_flatten(:,:,n) = ctf_correct_naive( all_intensity(:,:,n), defocus_vals(n), pixels, 1);
end
mean_intensity_ctf_correct_flatten = mean(all_intensity_ctf_correct_flatten(:,:,:),3);

all_intensity_ctf_correct_phase_flip = [];
for n = 1:size(all_intensity,3)
    all_intensity_ctf_correct_phase_flip(:,:,n) = ctf_correct_phase_flip( all_intensity(:,:,n), defocus_vals(n), pixels);
end
mean_intensity_ctf_correct_phase_flip = mean(all_intensity_ctf_correct_phase_flip(:,:,:),3);

%%
subplot(3,2,1);
show_map(imag(mean(all_amplitude(:,:,:),3)), pixels, 0,clim );
title( sprintf('Actual amplitude', size(all_intensity,3)) );


subplot(3,2,2);
show_map(mean_intensity_ctf_correct, pixels, 0,clim );
title( 'Intensity averaged,\newlinethen CTF-corrected naively with single defocus' );

subplot(3,2,3);
show_map(mean_intensity_ctf_correct_naive, pixels, 0,clim );
title( sprintf('Intensity CTF-corrected with naive', size(all_intensity,3)) );

subplot(3,2,4);
show_map( mean_intensity_ctf_correct_flatten, pixels, 0, clim);
title( sprintf('Intensity CTF-corrected with naive-flatten', size(all_intensity,3)) );

subplot(3,2,5);
show_map( mean_intensity_ctf_correct_phase_flip, pixels, 0, clim);
title( sprintf('Intensity CTF-corrected with phase-flip', size(all_intensity,3)) );

intensity_ctf_correct_weiner = ctf_correct_weiner( all_intensity, defocus_vals, pixels);
subplot(3,2,6);
show_map( intensity_ctf_correct_weiner, pixels, 0, clim);
title( sprintf('Intensity CTF-corrected with Weiner', size(all_intensity,3)) );

