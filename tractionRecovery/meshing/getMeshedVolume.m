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
function [node,elem0,elem,faces,mesh_options] = getMeshedVolume(mask,Lx,Ly,Lz,mesh_options,write_options)

%% Mask extension

%Extend the mask to have the same resolution in every direction

if Lz>Lx
    mask_extended = logical(imresize3(double(mask),[size(mask,1), size(mask,2), size(mask,3)*round(Lz/Lx)],'nearest'));
elseif Lz<Lx
    mask_extended = logical(imresize3(double(mask),[size(mask,1)*round(Lx/Lz),size(mask,2)*round(Lx/Lz),size(mask,3)],'nearest'));
else
    mask_extended = mask;
end

%% Processing and meshing

%Define the size for morphological operations
size_morphological_operator = 2;

%Pre-process the mask
mask_extended = preProcessMask(mask_extended,size_morphological_operator);

%Some random times it fails due to the iso2mesh functions so we keep
%trying until it works (with a maximum of 10 times)
error_count = 0;
iteration_count = 0;
while (iteration_count==error_count) && (iteration_count<=10)
    try
        [node,~,faces,elem0,elem] = meshVolume(mask_extended,mesh_options,Lx,Ly,Lz);
    catch
        error_count = error_count+1;
        disp(['Incompatible mesh: repeating (attempt ' num2str(error_count) ')...']);
    end
    iteration_count = iteration_count+1;
end
if error_count>10
    disp('THE MESH IS INCOMPATIBLE!')
end

%Correct the coordinates of the node according to the resolution
% min_res = min([Lx,Lx,Lz]);
% node = [min_res*node(:,1) min_res*node(:,2) min_res*node(:,3)];
if Lz>Lx
    node = [Lx*node(:,1) Ly*node(:,2) Lz/round(Lz/Lx,0)*node(:,3)];
elseif Lz<Lx
    node = [Lx/round(Lx/Lz,0)*node(:,1) Ly/round(Ly/Lz,0)*node(:,2) Lz*node(:,3)];
else
    node = [Lx*node(:,1) Ly*node(:,2) Lz*node(:,3)];
end

if not(isempty(write_options.write_path))  
    
    fid = fopen([write_options.write_path filesep write_options.write_filename '.flavia.msh'],'wt');
    fprintf(fid,'MESH  "Sprout" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 4');
    fprintf(fid,'\n','');
    fprintf(fid,'COORDINATES');
    fprintf(fid,'\n','');    
    
    rd=[(1:size(node,1))', node];
    fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
    
    fprintf(fid,'END COORDINATES');
    fprintf(fid,'\n','');
    fprintf(fid,'ELEMENTS');
    fprintf(fid,'\n','');       
    
    rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
    fprintf(fid,'%7i %7i %7i %7i %7i\n',rr');    
    
    fprintf(fid,'END ELEMENTS');
    fprintf(fid,'\n','');
    
    fprintf(fid,'MESH  "Gel" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 4');
    fprintf(fid,'\n','');
    fprintf(fid,'COORDINATES');
    fprintf(fid,'\n','');
    fprintf(fid,'END COORDINATES');
    fprintf(fid,'\n','');
    fprintf(fid,'ELEMENTS');
    fprintf(fid,'\n','');
    
    rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
    fprintf(fid,'%7i %7i %7i %7i %7i\n',rr');    
    
    fprintf(fid,'END ELEMENTS');
    fprintf(fid,'\n','');
    
    fclose(fid);
    
end

%% DEPRECATED

% mask_extended = false([size(mask,1),size(mask,2),size(mask,3)*round(Lz/Lx)]);
% for kk=1:size(mask,3)
%     for jj=1:round(Lz/Lx)
%         mask_extended(:,:,round(Lz/Lx)*(kk-1)+jj) = mask(:,:,kk);
%     end
% end
% 
% if Lz>Lx    
%     mask_extended = false([size(mask,1),size(mask,2),size(mask,3)*round(Lz/Lx)]);
%     for kk=1:size(mask,3)
%         for jj=1:round(Lz/Lx)
%             mask_extended(:,:,round(Lz/Lx)*(kk-1)+jj) = mask(:,:,kk);
%         end
%     end
%     
% elseif Lz<Lx
%     mask_extended = false([size(mask,1)*round(Lx/Lz),size(mask,2)*round(Lx/Lz),size(mask,3)]);
%     for ii=1:size(mask,1)
%         for jj=1:round(Lx/Lz)
%             mask_extended(round(Lx/Lz)*(ii-1)+jj,:,:) = mask(ii,:,:);
%         end
%     end
%     
%     for ii=1:size(mask,2)
%         for jj=1:round(Ly/Lz)
%             mask_extended(:,round(Lx/Lz)*(ii-1)+jj,:) = mask(:,ii,:);
%         end
%     end
% else
%     mask_extended = mask;
% end
