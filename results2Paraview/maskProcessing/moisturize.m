function [varargout] = moisturize (mask,res,filter_options)
% MOISTURIZE  Binary mask moisturizer (smoothing function).
%   MASK_SMOOTHED = MOISTURIZE(MASK,FILTER_OPTIONS,PLOT_FLAG) computes a
%   smoothed version MASK_SMOOTHED of the binary mask MASK.
%
%   MASK_SMOOTHED = MOISTURIZE(MASK,FILTER_OPTIONS) applies the filter
%   options available for function smoothpatch (check
%   https://nl.mathworks.com/matlabcentral/fileexchange/26710-smooth-triangulated-mesh).
%
%   MASK_SMOOTHED = MOISTURIZE(_,PLOT_FLAG) plots (plot_flag=1) the
%   interesting results.
%
%   [MASK_SMOOTHED,FV_SMOOTHED] = MOISTURIZE(...) returns both the smoothed
%   binary mask and the smoothed isosurface.
%
%   [_,FV_ORIGINAL] = MOISTURIZE(...) returns the isosurface of the
%   original binary mask.
%
% This function uses external functions SMOOTHPATCH and PLYGON2VOXEL from 
% Dirk-Jan Kroon (https://nl.mathworks.com/matlabcentral/profile/authors/1097878-dirk-jan-kroon)
% This function is written by Jorge Barrasa Fano (2018) KU Leuven.

% See also SMOOTHPATCH, ISOSURFACE, PLYGON2VOXEL PATCH
%% Check the inputs
if nargin == 2
    filter_options.mode = 1;
    filter_options.itt = 7;
    filter_options.lambda = 1;
    filter_options.sigma = 1;
    plot_flag = 0;
elseif nargin == 3
    plot_flag = 0;
end

%% Generate isosurface
field_size = size(mask);

[x,y,z] = meshgrid(1:field_size(2),1:field_size(1),1:field_size(3));
x = x*res(1);
y = y*res(2);
z = z*res(3);

fv_original = isosurface(x,y,z,mask,0.8);

%% Smooth the isosurface (Handle RAM bug)
pause(3);
try
    fv_smoothed = smoothpatch(fv_original,filter_options.mode,filter_options.itt,filter_options.lambda,filter_options.sigma);
catch exception
    fv_smoothed = smoothpatch(fv_original,filter_options.mode,filter_options.itt,filter_options.lambda,filter_options.sigma);
end

%% Convert to binary
mask_smoothed = polygon2voxel(fv_smoothed,size(mask),'none');
mask_smoothed = imfill(mask_smoothed,'holes'); %Fill all the volume inside the surface
%% Set output parameters
if nargout>=1
    varargout(1) = {mask_smoothed};
end
if nargout>=2
    varargout(2) = {fv_smoothed};
end
if nargout==3
    varargout(3) = {fv_original};
end
% %% Plot results
% if plot_flag
%     figure
%     subplot(131)
%     patch(fv_original,'facecolor',[255 127 127]/255,'edgecolor','none'), camlight;
%     daspect([1 1 1]),lighting flat
%     title('Original isosurface')
%     subplot(132)
%     patch(fv_smoothed,'facecolor',[255 127 127]/255,'edgecolor','none'), camlight;
%     daspect([1 1 1]),lighting flat
%     title('Smoothed isosurface')
%     %Convert the smoothed mask to check differences
%     fv_recovered = isosurface(mask_smoothed,0.8);
%     subplot(133)
%     patch(fv_recovered,'facecolor',[255 127 127]/255,'edgecolor','none'), camlight;
%     daspect([1 1 1]),lighting flat
%     title('Recovered isosurface from the smoothed mask')
% end