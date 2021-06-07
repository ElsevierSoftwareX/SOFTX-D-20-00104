%Run this file after stopping in function meshVolume after the line of v2s
sizes = [1.01 1.5 2 2.5 3 3.5];
container = [];
g=[];
for jj = 1:length(sizes)
    
    [mask_fv.vertices,mask_fv.faces,regions,holes]=v2s(mask_extended,0.5,sizes(jj),'cgalmesh');
        
    % Convert the points to um
    % points_um = points .* repmat([Lx Ly Lz/round(Lz/Lx,0)],[size(points,1),1]);
    points_um = mask_fv.vertices;
    %Get all pairs of adjacent points
    adjacent_points = unique([mask_fv.faces(:,[1,2]);...
        mask_fv.faces(:,[2,3]);...
        mask_fv.faces(:,[1,3])],'rows');
    
    %Compute the distances between adjacent points
    d = zeros(size(adjacent_points,1),1);
    for ii = 1:size(adjacent_points,1)
        d(ii) = pdist([points_um(adjacent_points(ii,1),:);points_um(adjacent_points(ii,2),:)]);
    end
    
    container=[container;d];
    g=[g;repmat({['s=' num2str(sizes(jj))]},size(d,1),size(d,2))];
end

figure
boxplot(container,g)

try
    mask_fv=smoothpatch(mask_fv,1,2);
catch exception
    mask_fv=smoothpatch(mask_fv,1,2);
end
 [node,el1,faces]=surf2mesh(mask_fv.vertices,mask_fv.faces,[-1 -1 -1],[size(mask_extended,1)+1 size(mask_extended,2)+1 size(mask_extended,3)+1],1,5000,regions,holes,1);

 %Get the elements corresponding to every region
elem0=el1((el1(:,5)==0),:);
elem=el1((el1(:,5)==1),:);

%Correct the coordinates of the node according to the resolution
node = [Lx*node(:,1) Ly*node(:,2) Lz/round(Lz/Lx,0)*node(:,3)];
 
figure
plotmesh(node,faces,el1);
title({'PARAMETERS :',...
    ['Gel meshing = ' num2str(mesh_options.maxvol_gel)],...
    ['Cell meshing = ' num2str(mesh_options.maxsize_cell)],...
    'RESULTANT ELEMENTS :',...
    ['N_{elements} hydrogel = ' num2str(size(elem0,1))]...
    ['N_{elements} cell = ' num2str(size(elem,1))],...
    ['N_{elements} TOTAL = ' num2str(size(elem,1)+size(elem0,1)) ] });