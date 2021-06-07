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

function results2paraview(path_names,experiment_name,folder_names,file_names,code_names)
f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Generating output files';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');

% Get the spatial resolution
res = getResolution(folder_names.tiff_raw);

%Define a pad size (necesary for exporting to Paraview files
pad_size = 3;
%% Get a ROI
% Get the channel names
folders = dir([path_names.store_path filesep experiment_name filesep folder_names.displacements]);
folders = folders(3:end);
tmp = false(1,length(folders));
for ii = 1:length(folders)
    if not(folders(ii).isdir)
        tmp(ii) = true;
    end
end
folders(tmp) = [];
clear tmp;

%Check if there is cell and nuclei (I only check it for the first folder (channel))
%since if there is cell and nuclei in one channel, there will also be for the other one
flag_cell = exist([path_names.store_path filesep experiment_name filesep folder_names.post_shift_corrected filesep folders(1).name filesep folder_names.cellSeg],'dir');

%Check for 3D
flag_3D = length(res) > 2;

%Make space in the RAM and crop the rest of the variables
save_path  = [path_names.store_path filesep experiment_name filesep folder_names.results_for_paraview];
if not(exist(save_path,'dir'))
    mkdir(save_path);
end


%Define the percentage of random points to sample from displacement field
random_percent = 0.2;

% Loop through the channels
for ii = 1:length(folders)
    %Get the name of the channel
    channel_name = folders(ii).name;
    
    %Get the time points
    files = dir([folders(ii).folder filesep channel_name]);
    files = files(3:end);    
    
    %Load the displacement field at the first time point
    load([folders(ii).folder filesep channel_name filesep files(1).name],'dispField');
    
    %Compute the magnitude
    if flag_3D
        mag_disp_field = sqrt(dispField.X.^2+dispField.Y.^2+dispField.Z.^2);
        %Max projection
        mag_disp_field_proj = max(mag_disp_field, [], 3);
    else
        mag_disp_field_proj = sqrt(dispField.X.^2+dispField.Y.^2);
    end
    clear dispField;
    %Plot
    figure;
    imshow(mag_disp_field_proj,[]);colormap(jet(256));colorbar;
    title(['Channel ' channel_name ': Please select a rectangular ROI and then double click in the rectangle']);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    
    %Load the cell mask if exists
    if flag_cell
        load([path_names.store_path filesep experiment_name filesep folder_names.post_shift_corrected filesep channel_name filesep folder_names.cellSeg filesep files(1).name]);
        if flag_3D
            cell_seg = max(cell_seg, [], 3);
        end
        hold on;
        visboundaries(cell_seg,'Color','w');
    end
    clear cell_seg;
    
    %Ask for a rectangle crop
    [~,crop_rect]=imcrop;
    close;
    
    %Crop the magnitude and add the padding just to get which are the final dimensions
    mag_disp_field = mag_disp_field(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
    mag_disp_field = padarray(mag_disp_field,pad_size*ones(1,length(res)),'replicate');
    %Get the field size
    field_size = size(mag_disp_field);
    clear mag_disp_field;
    
    %Get the number of random points that we need to sample
    nPts = round((random_percent/100)*prod(field_size));
    
    %Get the random indexes
    rand_indx = randperm(prod(field_size));
    sampling_indx = rand_indx(1:min(nPts,prod(field_size)));
    clear rand_indx nPts;
    
    if not(exist([save_path filesep channel_name],'dir'))
        mkdir([save_path filesep channel_name]);
    end
    
    %Create a container for the
    %Loop through the time points
    for tt = 1:length(files)
        %Update progress bar
        d.Value = (tt*ii)/(length(files)*length(folders));
        d.Message = sprintf(['Exporting files for ' channel_name ' at ' files(tt).name(1:end-4) '...']);
            
        %Load the displacement field
        load([files(tt).folder filesep files(tt).name],'dispField');        

        %Crop the displacement field
        dispField.X = dispField.X(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
        dispField.Y = dispField.Y(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
        if flag_3D
            dispField.Z = dispField.Z(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
        end
        
        %Add the pad
        dispField.X = padarray(dispField.X,pad_size*ones(1,length(res)),'replicate');
        dispField.Y = padarray(dispField.Y,pad_size*ones(1,length(res)),'replicate');
        if flag_3D
            dispField.Z = padarray(dispField.Z,pad_size*ones(1,length(res)),'replicate');
        end
        %% Generate .vtk files for Paraview     
        
        %Export the displacement filed
        disp2vtk(dispField,res,sampling_indx,channel_name,[save_path filesep channel_name filesep 'disp_' files(tt).name(1:end-4) '.vtk']);
        
        if flag_cell 
            
            %Load the cell mask
            load([path_names.store_path filesep experiment_name filesep folder_names.post_shift_corrected filesep channel_name filesep folder_names.cellSeg filesep files(tt).name]);
            %Crop the cell mask
            cell_seg = cell_seg(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
            %Add padding
            cell_seg = padarray(cell_seg,pad_size*ones(1,length(res)));
            %Export the cell mask
            cell2vtk(cell_seg,dispField,res,[save_path filesep channel_name filesep 'cell_' files(tt).name(1:end-4) '.vtk']);
            clear cell_seg;
            
        end
        
        %Nuclei
        if exist([path_names.store_path filesep experiment_name filesep folder_names.post_shift_corrected filesep channel_name filesep folder_names.nucleiSeg filesep files(tt).name],...
                'file')
            
            %Load the nuclei mask 
            load([path_names.store_path filesep experiment_name filesep folder_names.post_shift_corrected filesep channel_name filesep folder_names.nuclei filesep files(tt).name]);
            %Crop the cell mask
            nuclei_seg = nuclei_seg(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
            %Add padding
            nuclei_seg = padarray(nuclei_seg,pad_size*ones(1,length(res)));
            %Export the nuclei mask
            cell2vtk(nuclei_seg,dispField,res,[save_path filesep channel_name filesep 'nuclei_' files(tt).name(1:end-4) '.vtk']);
        
            clear nuclei_seg;
        end
        
        
        clear dispField;   
        
        
        %Load the Jacobian matrix (if it doesn't exist it will just throw a
        %warning)
        load([files(tt).folder filesep files(tt).name],'jacobianMat');
        
        if exist('jacobianMat','var')
            if not(isempty(jacobianMat))
                %Compute the volumetric strain and crop it
                vol_strain = tensorDeterminant(jacobianMat);
                clear jacobianMat;
                vol_strain = vol_strain(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
                
                %Add padding
                vol_strain = padarray(vol_strain,pad_size*ones(1,length(res)),'replicate');
                
                %Export the strain
                strain2vtk(vol_strain,res,sampling_indx,channel_name,[save_path filesep channel_name filesep 'strain_' files(tt).name(1:end-4) '.vtk']);
                clear vol_strain;
            end
        end
    end    
    
    save([save_path filesep channel_name filesep file_names.paraview_param],'sampling_indx','crop_rect','pad_size');
end
%Copy the template 
copyfile([code_names.results2Paraview filesep file_names.paraview_template],save_path);

% Close the dialog box
close(d);
delete(f);
