function show_map(intensity, pixels, dark_field, clim);
% show_map(intensity, pixels);
%

if ~exist( 'dark_field','var') dark_field = 0; end;
if ~exist( 'clim','var') 
    max_contrast= 0.2;
    clim = offset+max_contrast*[-1 1/2];
    if dark_field; clim = [0 0.1*max_contrast^2]; end;
end;


cla;


imagesc( pixels/1e-10,pixels/1e-10,intensity',clim);
%cmap = redwhiteblue(-max_contrast, max_contrast);
cmap = gray(256);
cmap = cmap(:,end:-1:1);
colormap(cmap);
colorbar();
title( 'Contrast, ÅxÅ pixels')
xlabel('x (Å)');ylabel('y (Å)')
axis image
