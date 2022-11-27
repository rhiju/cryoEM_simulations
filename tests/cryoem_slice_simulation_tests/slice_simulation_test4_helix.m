% 10 nm x 10 nm, 10 nm thickness
% pixel size 1 Å?
max_x = 50; max_z = 250; % in Angstroms
pixels   = [-max_x:1:max_x]*1e-10;
%pixels_z = [-50:1:50]*1e-10;
pixels_z = [-max_z:1:max_z]*1e-10;
[X,Y,Z] = ndgrid(pixels,pixels,pixels_z);

%
prefactor = (8.99e9)*(1.60217663e-19); % e/(4 pi epsilon_0) --> V*m
%
pdbstruct = pdbread( 'helix_10bp.pdb');
helix_grid = 0*X;
lambda_mu = 1e-10; % 1 Angstrom -- how well-screened the nucleus charge is.
for i = 1:length(pdbstruct.Model.Atom)
    % later update Z to be scattering factor
    switch pdbstruct.Model.Atom(i).element
        case 'H'
            Zq = 1;
        case 'C'
            Zq = 6;
        case 'N'
            Zq = 7;
        case 'O'
            Zq = 8;
        case 'P'
            Zq = 15;
        otherwise
            printf( 'Unknown! %s\n',pdbstruct.Model.Atom(i).element)
    end
    z = pdbstruct.Model.Atom(i).X;
    x = pdbstruct.Model.Atom(i).Y;
    y = pdbstruct.Model.Atom(i).Z;
    helix_grid( round(x+max_x)+1,...
                 round(y+max_x)+1,...
                 round(z+max_z)+1 ...
                 ) = Zq;
end
potential_helix_point = helix_grid * (4*pi) * prefactor * (lambda_mu)^2/(1e-10)^3; 
potential_helix = imgaussfilt3(potential_helix_point,1); % Use 1 Å smoothing

% delta functions for potentials of water.
h2o_per_A3 = 55 * 6.022e23 * 1e-27;
h2o_grid = (rand(size(X,1),size(X,2),size(X,3))>exp(-h2o_per_A3)); % this is a constant phase shift.
h2o_grid( find(potential_helix_point>0.5) ) = 0;
Z_h2o = 10; 
% summing over z later should give 4 * pi * lambda_mu^2 Z e /(4 pi epsilon0) , with
% units of V-m
potential_h2o = h2o_grid * Z_h2o * (4*pi) * prefactor * (lambda_mu)^2/(1e-10)^3; 
potential_h2o = imgaussfilt3(potential_h2o,1); % Use 1 Å smoothing


% potential for a dipole
potential_charge = 0*X;
n_phosphate = 0; n_counterion = 0;
for i = 1:length(pdbstruct.Model.Atom)
    % later update Z to be scattering factor
    if strcmp(pdbstruct.Model.Atom(i).AtomName ,'OP1') || ...
        strcmp(pdbstruct.Model.Atom(i).AtomName ,'OP2')
        Z_charge = -0.5; 
        dielectric = 1; %80;
        z = pdbstruct.Model.Atom(i).X;
        x = pdbstruct.Model.Atom(i).Y;
        y = pdbstruct.Model.Atom(i).Z;
        x1 = 1e-10*(x); 
        y1 = 1e-10*(y);
        z1 = 1e-10*(z);
        potential_charge = potential_charge + ...
            prefactor * (Z_charge/dielectric) * (  1./sqrt((X-x1).^2+(Y-y1).^2+(Z-z1).^2));
        n_phosphate = n_phosphate+Z_charge;
    end

    debye_length = 13e-10;
    if strcmp(pdbstruct.Model.Atom(i).AtomName ,'OP1') 
        Z_charge = 1; 
        dielectric = 1; %80;
        % This is a hack Gaussian -- would be much better to use
        % exponential, or even better to use NLPB solution to place
        % counterions.
        x2 = x1 + randn(1)*debye_length; 
        y2 = y1 + randn(1)*debye_length;
        z2 = z1 + randn(1)*debye_length;
        potential_charge = potential_charge + ...
            prefactor * (Z_charge/dielectric) * (  1./sqrt((X-x2).^2+(Y-y2).^2+(Z-z2).^2));
        n_counterion = n_counterion+Z_charge;
    end

end
fprintf( 'Total charge for phosphate, %f\n',n_phosphate);
fprintf( 'Total charge for counterion, %f\n',n_counterion);

% sum of potentials
potential =  potential_helix + potential_charge + potential_h2o - h2o_per_A3* Z_h2o * (4*pi) * prefactor * (lambda_mu)^2/(1e-10)^3;
proj_potential = sum(potential,3)*1e-10;
max(max(proj_potential))

%
% Phase shift
E = 300e3; % in eV
lambda = 2.24e-12; % 300 keV electron -->2.24 pm 
phase_shift = pi * proj_potential/lambda/E;
% fprintf( 'Maximum phase shift in 1 Å pixel is: %f rad\n',max(max(phase_shift)));
% fprintf( 'Integrated phase shift in 1 Å pixel is: %f rad * Å^2\n',sum(sum(phase_shift)));

contrast = 2 * phase_shift;

% show it!
figure(1)
max_contrast= 10;
imagesc( pixels,pixels,contrast',max_contrast*[-1 1]);
cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
title( 'Contrast, ÅxÅ pixels')
axis image

figure(2)
max_contrast= 2;
contrast_smooth = imgaussfilt(contrast,5);
imagesc( pixels,pixels,contrast_smooth',max_contrast*[-1 1]);
colormap(cmap);
colorbar();
title( 'Contrast, smoothed, ÅxÅ pixels')
axis image


%
figure(3)
clf
potential_show =  potential_helix + potential_charge;
p = patch(isosurface(pixels,pixels,pixels_z, potential_show, -6));
p.FaceColor = 'red';
p.EdgeColor = 'none';
p.FaceAlpha = 1;
p = patch(isosurface(pixels,pixels,pixels_z, potential_show, 6));
p.FaceColor = 'blue';
p.EdgeColor = 'none';
p.FaceAlpha = 1;
daspect([1 1 1])
view(3)
camlight; lighting phong
