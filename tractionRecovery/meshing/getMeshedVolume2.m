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

function [node,elem0,elem,faces,mesh_options] = getMeshedVolume2(mask,Lx,Ly,Lz,mesh_options,mech_props,write_options,dispField)
%% USING DISPLACEMENT FIELD

%Extend the mask in case the resolutions are different
if Lz>Lx
    mask_extended = logical(imresize3(double(mask),[size(mask,1), size(mask,2), size(mask,3)*round(Lz/Lx)],'nearest'));
elseif Lz<Lx
    mask_extended = logical(imresize3(double(mask),[size(mask,1)*round(Lx/Lz),size(mask,2)*round(Lx/Lz),size(mask,3)],'nearest'));
else
    mask_extended = mask;
end
size_morphological_operator = 2;
mask_extended = preProcessMask(mask_extended,size_morphological_operator);


%Some random times it fails due to the iso2mesh functions so we keep
%trying until it works (with a maximum of 10 times)
error_count = 0;
iteration_count = 0;
while (iteration_count==error_count) && (iteration_count<=10)
    try
        %Do the FE mesh
        [mask_fv.vertices,mask_fv.faces,regions,holes] = v2s(mask_extended,0.5,mesh_options.maxsize_cell,'cgalmesh');
           
        %Make the rows be Y and the columns be X
        mask_fv.vertices = [mask_fv.vertices(:,2) mask_fv.vertices(:,1) mask_fv.vertices(:,3)];
        regions = [regions(:,2) regions(:,1) regions(:,3)];
        if not(isempty(holes))
            holes = [holes(:,2) holes(:,1) holes(:,3)];
        end
        [mask_fv.vertices,mask_fv.faces] = meshcheckrepair(mask_fv.vertices,mask_fv.faces);
        
        %Smooth the surface
        if mesh_options.post_smoothing >0
            conn=meshconn(mask_fv.faces(:,1:3),size(mask_fv.vertices,1));
            mask_fv.vertices=smoothsurf(mask_fv.vertices,[],conn,mesh_options.post_smoothing,mesh_options.post_smoothing_amount,mesh_options.post_smoothing_method);
            [mask_fv.vertices,mask_fv.faces]=meshcheckrepair(mask_fv.vertices,mask_fv.faces);
        end        
        
        %Resample
%         [mask_fv.vertices,mask_fv.faces]=meshresample(mask_fv.vertices,mask_fv.faces(:,1:3),0.95);
        
        %Find the maximum of the magnitude of the displacement field in the surface      
        mag_disp = sqrt(dispField.X.^2+dispField.Y.^2+dispField.Z.^2);
        if Lz>Lx
            mag_disp_extended = (imresize3(double(mag_disp),[size(mag_disp,1), size(mag_disp,2), size(mag_disp,3)*round(Lz/Lx)],'cubic'));
        elseif Lz<Lx
            mag_disp_extended = (imresize3(double(mag_disp),[size(mag_disp,1)*round(Lx/Lz),size(mag_disp,2)*round(Lx/Lz),size(mag_disp,3)],'cubic'));
        else
            mag_disp_extended = mag_disp;
        end
        [Gmag,~,~] = imgradient3(mag_disp_extended); 
        thresh = multithresh(Gmag,4);
        Gmag_seg = imquantize(Gmag,thresh);
        Gmag_seg2 = zeros(size(Gmag_seg));
        Gmag_seg2(Gmag_seg==1) = 8;
        Gmag_seg2(Gmag_seg==2) = 4;
        Gmag_seg2(Gmag_seg==3) = 2;
        Gmag_seg2(Gmag_seg==4) = 1;
        Gmag_seg2(Gmag_seg==5) = 0.1;
        clear Gmag_seg;
        
        %Define maximum tetrahedra element volume
%         mesh_options.quality = 1.414;
        minangle = 18;
        maxangle = 100;
%         mesh_options.maxvol_gel = 200;
        ISO2MESH_TETGENOPT = ['-A -pq' num2str(mesh_options.quality) '/' num2str(minangle) ...
                      ' -a' num2str(mesh_options.maxvol_gel) ' -V -o/' num2str(maxangle)]; 
        [node,el1,faces]=surf2mesh(mask_fv.vertices,mask_fv.faces,[-1 -1 -1],[size(mask_extended,2)+1 size(mask_extended,1)+1 size(mask_extended,3)+1],1,mesh_options.maxvol_gel,regions,holes,1,'tetgen1.5');
        %Change the label of the gel elements to 0 instead of 2
        el1(not(el1(:,5)==1),5) = 0;
        
        
        if mesh_options.refine_flag
            %Define interpolants
            field_size = size(Gmag);
            [x,y,z] = meshgrid((1:field_size(2))-1,(1:field_size(1))-1,(1:field_size(3))-1);
            F = griddedInterpolant(y,x,z,Gmag_seg2,'nearest'); % Build the interpolation function
            Gmag_nodes = F(node(:,2),node(:,1),node(:,3));
            %         Gmag_nodes(Gmag_nodes<0) = 0;
            %Normalize the gradient to the specified range
            %         range = [0.2 2];
            %         edge_sizes = rescale(-Gmag_nodes,0.2,2);
            edge_sizes = Gmag_nodes;
            %[newnode,newelem,newface]=meshrefine(node,el1,faces,edge_sizes);
            [node,el1,faces]=meshrefine(node,el1,faces,edge_sizes);
        end
%         figure
%         plotmesh(newnode,newface);
%         figure
%         plotmesh(node,faces)
% 
%         figure
%         hsurf=trimesh(faces(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');
%         hsurf.EdgeColor = 'none';
%         hsurf.EdgeAlpha = 0.2;
%         hold on;
%         %         z0=cpoint(3);
%         z0=mean(node(:,3));
%         z0 = 20;
%         plane=[min(node(:,1)) min(node(:,2)) z0
%             min(node(:,1)) max(node(:,2)) z0
%             max(node(:,1)) min(node(:,2)) z0];
%         [cutpos,cutvalue,facedata]=qmeshcut(el1(:,1:4),node(:,1:3),node(:,1),plane);
%         hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','g');
% %         scatter3(cpoint(1),cpoint(2),cpoint(3),80,'r','filled')
%         daspect([1,1,1])
%         title(['#Elements = ' num2str(size(el1,1)) '; #Nodes = ' num2str(size(node,1))...
%             '; Critical point size = ' num2str(mesh_options.size_critical_points) '; Cell surface size = ' num2str(mesh_options.maxsize_cell)]);

%         
%         figure
%         newnode,newelem,newface
%         hsurf=trimesh(newface(:,1:3),newnode(:,1),newnode(:,2),newnode(:,3),'facecolor','none');
%         hsurf.EdgeColor = 'w';
%         hsurf.EdgeAlpha = 1;
%         hold on;
%         z0=mean(newnode(:,3));
%         plane=[min(newnode(:,1)) min(newnode(:,2)) z0
%             min(newnode(:,1)) max(newnode(:,2)) z0
%             max(newnode(:,1)) min(newnode(:,2)) z0];
%         [cutpos,cutvalue,facedata]=qmeshcut(newelem(:,1:4),newnode(:,1:3),newnode(:,1),plane);
%         hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','g');
%         daspect([1,1,1])
%         title(['#Elements = ' num2str(size(el1,1)) '; #Nodes = ' num2str(size(node,1))...
%             '; Critical point size = ' num2str(mesh_options.size_critical_points) '; Cell surface size = ' num2str(mesh_options.maxsize_cell)]);

        
        %Correct the coordinates of the node according to the resolution
%         min_res = min([Lx,Lx,Lz]);
%         node = [min_res*node(:,1) min_res*node(:,2) min_res*node(:,3)];
%         cpoint = min_res *cpoint;
        if Lz>Lx
            node = [Lx*node(:,1) Ly*node(:,2) Lz/round(Lz/Lx,0)*node(:,3)];
%             cpoint = cpoint .* [Lx Ly Lz/round(Lz/Lx,0)];
        elseif Lz<Lx
            node = [Lx/round(Lx/Lz,0)*node(:,1) Ly/round(Ly/Lz,0)*node(:,2) Lz*node(:,3)];
%             cpoint = cpoint .* [Lx/round(Lx/Lz,0) Ly/round(Ly/Lz,0) Lz];
        else
            node = [Lx*node(:,1) Ly*node(:,2) Lz*node(:,3)];
%             cpoint = cpoint .* [Lx Ly Lz];
        end
        switch mech_props.abaqus.element_type
            case 'C3D4'
                %Get the elements corresponding to every region
                elem0=el1((el1(:,5)==0),:);
                elem=el1((el1(:,5)==1),:);
%                 q = meshquality(node,[elem0;elem]);
%                 save([write_options.write_path filesep write_options.write_filename '_qels.mat'],'q');
                
            case 'C3D10'
                [node,el1]= getC3D10elements(node,el1);
                %Get the elements corresponding to every region
                elem0=el1((el1(:,11)==0),:);
                elem=el1((el1(:,11)==1),:);
        end
        clear el1;
        
        if not(isempty(write_options.write_path))
            
            fid = fopen([write_options.write_path filesep write_options.write_filename '.flavia.msh'],'wt');
            switch mech_props.abaqus.element_type
                case 'C3D4'
                    fprintf(fid,'MESH  "Sprout" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 4');
                case 'C3D10'
                    fprintf(fid,'MESH  "Sprout" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 10');
            end
            fprintf(fid,'\n','');
            fprintf(fid,'COORDINATES');
            fprintf(fid,'\n','');
            
            rd=[(1:size(node,1))', node]; 
            
            switch mech_props.abaqus.element_type
                case 'C3D4'
                    fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
                case 'C3D10'
                    fprintf(fid,'%7i %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n',rd');
            end
            
            fprintf(fid,'END COORDINATES');
            fprintf(fid,'\n','');
            fprintf(fid,'ELEMENTS');
            fprintf(fid,'\n','');
            
            switch mech_props.abaqus.element_type
                case 'C3D4'
                    rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
                    fprintf(fid,'%7i %7i %7i %7i %7i\n',rr');
                case 'C3D10'
                    rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4),...
                        elem(:,5), elem(:,6), elem(:,7), elem(:,8), elem(:,9), elem(:,10)];
                    fprintf(fid,'%7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i \n',rr');
            end    
            
            fprintf(fid,'END ELEMENTS');
            fprintf(fid,'\n','');
            
            switch mech_props.abaqus.element_type
                case 'C3D4'
                    fprintf(fid,'MESH  "Gel" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 4');
                case 'C3D10'
                    fprintf(fid,'MESH  "Gel" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 10');
            end
            fprintf(fid,'\n','');
            fprintf(fid,'COORDINATES');
            fprintf(fid,'\n','');
            fprintf(fid,'END COORDINATES');
            fprintf(fid,'\n','');
            fprintf(fid,'ELEMENTS');
            fprintf(fid,'\n','');
            
            switch mech_props.abaqus.element_type
                case 'C3D4'
                    rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
                    fprintf(fid,'%7i %7i %7i %7i %7i\n',rr');                    
                case 'C3D10'
                    rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4),...
                        elem0(:,5), elem0(:,6), elem0(:,7), elem0(:,8), elem0(:,9), elem0(:,10)];
                    fprintf(fid,'%7i %7i %7i %7i %7i %7i %7i %7i %7i %7i %7i \n',rr');
            end
                        
            fprintf(fid,'END ELEMENTS');
            fprintf(fid,'\n','');
            
            fclose(fid);
            
        end
        
        
        
    catch
        error_count = error_count+1;
        disp(['Incompatible mesh: repeating (attempt ' num2str(error_count) ')...']);
    end
    iteration_count = iteration_count+1;
end
if error_count>10
    disp('THE MESH IS INCOMPATIBLE!')
end


% 
% 
% plane=[cpoint(1) min(node(:,2)) min(node(:,3))
%        cpoint(1) max(node(:,2)) min(node(:,3))
%        cpoint(1) min(node(:,2)) max(node(:,3))];
% [cutpos,cutvalue,facedata]=qmeshcut(el1(:,1:4),node(:,1:3),node(:,1),plane);
% hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','g');
% plane=[min(node(:,1)) cpoint(2) min(node(:,3))
%        min(node(:,1)) cpoint(2) max(node(:,3))
%        max(node(:,1)) cpoint(2) max(node(:,3))];
% [cutpos,cutvalue,facedata]=qmeshcut(el1(:,1:4),node(:,1:3),node(:,1),plane);
% hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','g');

%% Dilated version
% se = strel('sphere',10);
% mask_dil = imdilate(mask_extended,se);
% 
% figure
% h1 = PATCH_3Darray(mask_dil1);
% h1.FaceAlpha = 0.1;hold on;
% h2 = PATCH_3Darray(mask_extended);
% h2.FaceColor = 'r';
% 
% % Mix both masks
% mask_both = mask_dil+(mask_extended>0)+1;
% 
% opt(1).radbound=1; % brain surface element size bound
% opt(2).radbound=200; % head surface element size bound
% opt(3).radbound=200; % brain surface element size bound
% 
% [node2,elem2,face2]=vol2mesh(uint8(mask_both),1:size(mask_extended,1),1:size(mask_extended,2),1:size(mask_extended,3),opt,'0=20:1=20:2=20',1,'cgalmesh');
% size(elem2,1)
% figure
% plotmesh(node2(:,[2 1 3]),face2);
% 
% 
% figure;
% hsurf=trimesh(faces(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');
% hsurf.EdgeColor = 'k';
% hsurf.EdgeAlpha = 0.2;
% hold on;
% if(isoctavemesh)
%   hcut=patch('Faces',facedata,'Vertices',cutpos);
% else
%   hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','interp');
% end
% hcut.FaceColor='g';
% %set(hcut, 'linestyle','none')
% axis equal;
% 
% 
% 
% 
% z0=mean(node(:,3));
% 
% plane=[min(node(:,1)) min(node(:,2)) z0
%        min(node(:,1)) max(node(:,2)) z0
%        max(node(:,1)) min(node(:,2)) z0];
% 
% % run qmeshcut to get the cross-section information at z=mean(node(:,1))
% % use the x-coordinates as the nodal values
% 
% [cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node(:,1:3),node(:,1),plane);
% 
% figure;
% hsurf=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');
% hsurf.EdgeColor = 'k';
% hsurf.EdgeAlpha = 0.2;
% hold on;
% if(isoctavemesh)
%   hcut=patch('Faces',facedata,'Vertices',cutpos);
% else
%   hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','interp');
% end
% hcut.FaceColor='g';
% %set(hcut, 'linestyle','none')
% axis equal;
% 
% 
% 
% 
% [mask_fv.vertices,mask_fv.faces,regions,holes]=v2s(mask_extended,0.5,mesh_options.maxsize_cell,'cgalmesh');
% 
% 
