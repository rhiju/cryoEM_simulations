function counts = simulate_counts( intensity, dose, pixels)
% counts = simulate_counts( intensity, dose, pixels)
%
% Inputs:
%  intensity = 2D map of intensity (pass-through intensity = 1)
%  dose      = fluence in e-/Å^2
%  pixels    = pixel locations (1D) in m.
%
% Output:
%  counts = 2D map of counts (will be integers), Poisson-distributed
%
% (C) R. Das, Stanford (2022)

counts = [];
pixel_size = (pixels(2)-pixels(1))/1e-10;
pixel_area = pixel_size^2;
counts = poissrnd( intensity*dose*pixel_area );