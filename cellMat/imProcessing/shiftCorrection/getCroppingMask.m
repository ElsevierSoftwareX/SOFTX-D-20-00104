function cuboid_shape = getCroppingMask(field_size,axang,shift)
%Create an initial black mask
mask = ones(field_size);

%Get the angle (it is always in the 4th position, even in 2D)
D = rad2deg(axang(4));

%Transform the mask
if length(field_size)==2
    mask_tx = imrotate(mask,D,axang(1:2),'linear','crop');
else
    mask_tx = imrotate3(mask,D,axang(1:3),'linear','crop','FillValues',0);
end
mask_tx = imtranslate(mask_tx,shift);

% Convert to logical
mask_tx = mask_tx==1;
mask = logical(mask);

%% Create a central cube
%Create a padding before doin
n = 3;
% % Compute the distance transform
% dt_mask_tx = bwdist(not(padarray(mask_tx,[n n n])));
% %Remove the padding
% dt_mask_tx = dt_mask_tx(n+1:end-n,n+1:end-n,n+1:end-n);
% 
% % Find the max of the dt:
% max_distance = max(dt_mask_tx(:));

%Find centroid of the mask
centroid = regionprops3(mask_tx,'Centroid');
centroid = round(table2array(centroid));
%Create a container for a cuboid that will start expanding
cuboid_shape = false(size(mask));
cuboid_shape(centroid(2),centroid(1),centroid(3)) = true;

% tmp = dt_mask_tx == max_distance;
% cuboid_shape(tmp) = true;
% clear dt_mask_tx max_distance;

%Dilate the cuboid until the overlap is the same
se = strel('cube',3);
while sum(and(cuboid_shape(:),mask_tx(:))) == sum(cuboid_shape(:))
%     n_overlap_previous_it = n_overlap;
    cuboid_shape = imdilate(cuboid_shape,se);
%     dip_image(and(cuboid_shape,mask_tx))
%     n_overlap = sum(and(cuboid_shape(:),mask_rot(:)));
end
%Correct the cuboid in case we have gone out of the mask
cuboid_shape = imerode(cuboid_shape,se);
cuboid_shape = and(cuboid_shape,mask_tx);
% joinchannels('RGB',cuboid_shape,mask_tx)
clear se;
%% Extra extensions

%Get the boundaries of the object
n = 3;
cuboid_shape = padarray(cuboid_shape,[n n n]);
cuboid_boundaries = cuboid_shape & not(imerode(cuboid_shape,strel('cube',3)));
cuboid_shape = cuboid_shape(n+1:end-n,n+1:end-n,n+1:end-n);
cuboid_boundaries = cuboid_boundaries(n+1:end-n,n+1:end-n,n+1:end-n);
clear n;
%Find the coordinates
[rows,columns,slices] = ind2sub(size(cuboid_boundaries),find(cuboid_boundaries));
clear cuboid_boundaries;
%Get the 2 most common rows
mc_rows = mode(rows);
rows(rows==mc_rows) = [];
mc_rows(2) = mode(rows);
clear rows;

%Get the 2 most common columns
mc_columns = mode(columns);
columns(columns==mc_columns) = [];
mc_columns(2) = mode(columns);
clear columns;

%Get the 2 most common slices
mc_slices = mode(slices);
slices(slices==mc_slices) = [];
mc_slices(2) = mode(slices);
clear slices;

%Extend the two sides of each dimension
cuboid_shape = extendCube(cuboid_shape,mask_tx,mc_rows(1),'rows');
cuboid_shape = extendCube(cuboid_shape,mask_tx,mc_rows(2),'rows');
% joinchannels('RGB',cuboid_shape,mask_tx)
cuboid_shape = extendCube(cuboid_shape,mask_tx,mc_columns(1),'columns');
cuboid_shape = extendCube(cuboid_shape,mask_tx,mc_columns(2),'columns');
% joinchannels('RGB',cuboid_shape,mask_tx)
cuboid_shape = extendCube(cuboid_shape,mask_tx,mc_slices(1),'slices');
cuboid_shape = extendCube(cuboid_shape,mask_tx,mc_slices(2),'slices');
% joinchannels('RGB',cuboid_shape,mask_tx)

% %Crop the image following the mask
% im_rot = imtranslate(im,[10 0 0]);
% im_rot = imrotate3(im_rot,5,[0 0 1],'linear','crop','FillValues',255);
% 
% [rows, columns, slices] = ind2sub(size(cuboid_shape),find(cuboid_shape));
% tmp = im_rot(rows(1):rows(end), columns(1):columns(end), slices(1):slices(end));

end

%% Extra function

function cuboid_shape = extendCube(cuboid_shape,mask_rot,mc,direction)
cuboid_shape_dir = false(size(cuboid_shape));
switch direction
    case 'rows'
        cuboid_shape_dir(mc,:,:) = true;
        se_shape = [3 1 1];
    case 'columns'
        cuboid_shape_dir(:,mc,:) = true;
        se_shape = [1 3 1];
    case 'slices'
        cuboid_shape_dir(:,:,mc) = true;
        se_shape = [1 1 3];
end
cuboid_shape_dir = and(cuboid_shape_dir,cuboid_shape);

%Dilate the cuboid until the overlap is the same
se = strel('cuboid',se_shape);
% n_overlap = sum(and(cuboid_shape_dir(:),mask_rot(:)));
% n_overlap_previous_it = 0;
%
while (sum(and(cuboid_shape_dir(:),mask_rot(:))) == sum(cuboid_shape_dir(:))) && isequal(sum(imdilate(cuboid_shape_dir,se),'all'),sum(imdilate(cuboid_shape_dir,se,'full'),'all'))
%     n_overlap_previous_it = n_overlap;
    cuboid_shape_dir = imdilate(cuboid_shape_dir,se);
%     dip_image(and(cuboid_shape_dir,mask_rot))
%     n_overlap = sum(and(cuboid_shape_dir(:),mask_rot(:)));
end
% close all;
cuboid_shape_dir = imerode(cuboid_shape_dir,se);
cuboid_shape = and(or(cuboid_shape_dir,cuboid_shape),mask_rot);
end