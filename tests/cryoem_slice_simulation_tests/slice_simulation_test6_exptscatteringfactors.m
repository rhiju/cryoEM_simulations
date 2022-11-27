% 10 nm x 10 nm, 20 nm thickness
% pixel size 1 Å?
max_x = 60; max_z = 100; % in Angstroms
pixels   = [-max_x:1:max_x]*1e-10;
pixels_z = [-max_z:1:max_z]*1e-10;
[X,Y,Z] = ndgrid(pixels,pixels,pixels_z);

%
prefactor = (8.99e9)*(1.60217663e-19); % e/(4 pi epsilon_0) --> V*m
%
pdbstruct = pdbread( 'helix_40bp.pdb');
% Let's apply an offset
for i = 1:length(pdbstruct.Model.Atom)
    all_r(:,i) = [pdbstruct.Model.Atom(i).X,pdbstruct.Model.Atom(i).Y,pdbstruct.Model.Atom(i).Z];
end
mean_r = mean(all_r');
mean_x = mean_r(1);
mean_y = mean_r(2);
mean_z = mean_r(3);

% Prepare grid of scattering potentials from helix.
helix_grid = 0*X;
lambda_mu = 1e-10; % 1 Angstrom -- how well-screened the nucleus charge is.
for i = 1:length(pdbstruct.Model.Atom)
    % later update Z to be scattering factor
    switch pdbstruct.Model.Atom(i).element
        case 'H'
            A = 0.5e-10; 
        case 'C'
            A = 2.4e-10;
        case 'N'
            A = 2.2e-10;
        case 'O'
            A = 2.0e-10;
        case 'P'
            A = 5.5e-10;
        otherwise
            printf( 'Unknown! %s\n',pdbstruct.Model.Atom(i).element)
    end
    r = [pdbstruct.Model.Atom(i).X,pdbstruct.Model.Atom(i).Y,pdbstruct.Model.Atom(i).Z]-mean_r;
    r = circshift(r,1); % rotate 90°
    helix_grid( round(r(1)+max_x)+1,...
                 round(r(2)+max_x)+1,...
                 round(r(3)+max_z)+1 ...
                 ) = A;
end

E = 300e3; % in eV
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 
amplitude_to_volt_m3 = (E*lambda^2/pi);  % converts from scattering amplitude (in m) to Volt-m^3
potential_helix_point = helix_grid * amplitude_to_volt_m3/(1e-10)^3; 
potential_helix = imgaussfilt3(potential_helix_point,1); % Use 1 Å smoothing

% delta functions for potentials of water.
h2o_per_A3 = 55 * 6.022e23 * 1e-27;
%h2o_grid = (rand(size(X,1),size(X,2),size(X,3))>exp(-h2o_per_A3)); % this is a constant phase shift.
h2o_grid = (rand(size(X,1),size(X,2),size(X,3))<h2o_per_A3); % this is a constant phase shift.
h2o_grid( find(helix_grid>0.5) ) = 0;
Z_h2o = 10; 
% summing over z later should give 4 * pi * lambda_mu^2 Z e /(4 pi epsilon0) , with
% units of V-m
potential_h2o = h2o_grid * 3.0e-10*amplitude_to_volt_m3/(1e-10)^3; 
potential_h2o = imgaussfilt3(potential_h2o,1); % Use 1 Å smoothing


% potential for a dipole
potential_charge = 0*X;
n_phosphate = 0; n_counterion = 0;
Rcut = 3e-10; % Cutoff of full strength Coulomb. Outside this radius, assume water screens Coulomb potential to 0
for i = 1:length(pdbstruct.Model.Atom)
    % later update Z to be scattering factor
    if strcmp(pdbstruct.Model.Atom(i).AtomName ,'OP1') || ...
        strcmp(pdbstruct.Model.Atom(i).AtomName ,'OP2')
        Z_charge = 1*-0.5; 
        r = [pdbstruct.Model.Atom(i).X,pdbstruct.Model.Atom(i).Y,pdbstruct.Model.Atom(i).Z]-mean_r;
        r = circshift(r,1); % rotate 90°
        x1 = 1e-10*r(1);
        y1 = 1e-10*r(2);
        z1 = 1e-10*r(3);
        R = sqrt((X-x1).^2+(Y-y1).^2+(Z-z1).^2);
        potential_charge = potential_charge + ...
            prefactor * (Z_charge) * ((R<Rcut)./R);
        n_phosphate = n_phosphate+Z_charge;
    end

    debye_length = 13e-10;
    if strcmp(pdbstruct.Model.Atom(i).AtomName ,'OP1') 
        Z_charge = 1; 
        % This is a hack Gaussian -- would be much better to use
        % exponential, or even better to use NLPB solution to place
        % counterions.
        x2 = x1 + randn(1)*debye_length; 
        y2 = y1 + randn(1)*debye_length;
        z2 = z1 + randn(1)*debye_length;
        R = sqrt((X-x2).^2+(Y-y2).^2+(Z-z2).^2);
        potential_charge = potential_charge + ...
            prefactor * (Z_charge) * ((R<Rcut)./R);
        n_counterion = n_counterion+Z_charge;
    end

end
fprintf( 'Total charge for phosphate, %f\n',n_phosphate);
fprintf( 'Total charge for counterion, %f\n',n_counterion);

% sum of potentials
%potential =  potential_helix + potential_charge;
potential =  potential_helix + potential_charge + potential_h2o - h2o_per_A3* 3.0e-10*amplitude_to_volt_m3/(1e-10)^3;
%potential =  potential_helix + potential_h2o - h2o_per_A3* 3.0e-10*amplitude_to_volt_m3/(1e-10)^3;
proj_potential = sum(potential,3)*1e-10;
max(max(proj_potential))

%
% Phase shift
phase_shift = pi * proj_potential/lambda/E;
% fprintf( 'Maximum phase shift in 1 Å pixel is: %f rad\n',max(max(phase_shift)));
% fprintf( 'Integrated phase shift in 1 Å pixel is: %f rad * Å^2\n',sum(sum(phase_shift)));

contrast = 2 * phase_shift;

%
% show it!
subplot(2,2,1);
max_contrast= 0.2;
imagesc( pixels,pixels,contrast',max_contrast*[-1 1]);
cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
title( 'Contrast, ÅxÅ pixels')
axis image

subplot(2,2,2);
max_contrast= 0.2;
contrast_smooth = imgaussfilt(contrast,5);
imagesc( pixels,pixels,contrast_smooth',max_contrast*[-1 1]);
colormap(cmap);
colorbar();
title( 'Contrast, smoothed, ÅxÅ pixels')
axis image


%
subplot(2,1,2);
cla
potential_show =  potential_helix + potential_charge;
max_V = 2;
p = patch(isosurface(pixels,pixels,pixels_z, potential_show, -max_V));
p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = 1;
p = patch(isosurface(pixels,pixels,pixels_z, potential_show, max_V));
p.FaceColor = 'blue';
p.EdgeColor = 'none';
p.FaceAlpha = 1;
daspect([1 1 1])
view(3)
camlight; lighting phong
title( sprintf('Contours at +/- %5.1f V',max_V) );
