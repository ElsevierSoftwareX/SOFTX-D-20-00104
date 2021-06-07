%    This file is part of TFMLAB.
%     Copyright (C) 2020 MAtrix
%     Copyright (C) 2020 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2020 Esc. Téc. Sup. de Ingeniería, Universidad de Sevilla (Spain)
%
%     This library is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published
%     by the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This software is provided "as is",
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details
%     <http://www.gnu.org/licenses/>.

function [node,el1,faces,elem0,elem] = meshVolume(mask_extended,mesh_options,Lx,Ly,Lz)
%Transform the mask to tethrahedral elements
[mask_fv.vertices,mask_fv.faces,regions,holes]=v2s(mask_extended,0.5,mesh_options.maxsize_cell,'cgalmesh');

%Make the rows be Y and the columns be X
mask_fv.vertices = [mask_fv.vertices(:,2) mask_fv.vertices(:,1) mask_fv.vertices(:,3)];
regions = [regions(:,2) regions(:,1) regions(:,3)];
if not(isempty(holes))
    holes = [holes(:,2) holes(:,1) holes(:,3)];
end

%Check the average distance between vertices
% mean_dist = calculateMeanNodeDist(mask_fv.vertices,Lx,Ly,Lz,mask_fv.faces);
% 
% figure;
% plotmesh(mask_fv.vertices,mask_fv.faces);

%Smoothen the surface
% if mesh_options.post_smoothing>0
%     try
%         mask_fv=smoothpatch(mask_fv,1,mesh_options.post_smoothing);
%     catch exception
%         mask_fv=smoothpatch(mask_fv,1,mesh_options.post_smoothing);
%     end
% end
if mesh_options.post_smoothing >0
    conn=meshconn(mask_fv.faces(:,1:3),size(mask_fv.vertices,1));
    mask_fv.vertices=smoothsurf(mask_fv.vertices,[],conn,mesh_options.post_smoothing,mesh_options.post_smoothing_amount,mesh_options.post_smoothing_method);
    [mask_fv.vertices,mask_fv.faces]=meshcheckrepair(mask_fv.vertices,mask_fv.faces);
end

%Define maximum tetrahedra element volume
% [node,el1,faces]=surf2mesh(mask_fv.vertices,mask_fv.faces,[-1 -1 -1],[size(mask_extended,1)+1 size(mask_extended,2)+1 size(mask_extended,3)+1],1,mesh_options.maxvol_gel,regions,holes,1);
[node,el1,faces]=surf2mesh(mask_fv.vertices,mask_fv.faces,[-1 -1 -1],[size(mask_extended,2)+1 size(mask_extended,1)+1 size(mask_extended,3)+1],1,mesh_options.maxvol_gel,regions,holes,1);



%Get the elements corresponding to every region
elem0=el1((el1(:,5)==0),:);
elem=el1((el1(:,5)==1),:);


end

function mean_dist = calculateMeanNodeDist(points,Lx,Ly,Lz,faces)


% Convert the points to um
points_um = points .* repmat([Lx Ly Lz/round(Lz/Lx,0)],[size(points,1),1]);

%Get all pairs of adjacent points
adjacent_points = unique([faces(:,[1,2]);...
    faces(:,[2,3]);...
    faces(:,[1,3])],'rows');

%Compute the distances between adjacent points (this approach is faster
%than cellfun!) 
d = zeros(size(adjacent_points,1),1);
for ii = 1:size(adjacent_points,1)
    d(ii) = pdist([points_um(adjacent_points(ii,1),:);points_um(adjacent_points(ii,2),:)]);
end

%Compute the minimum distance
mean_dist = mean(d);

end
