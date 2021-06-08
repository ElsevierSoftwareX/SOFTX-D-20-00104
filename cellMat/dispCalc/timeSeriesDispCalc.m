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

function timeSeriesDispCalc(path_names,experiment_name,folder_names,dispParam,disp_param_file_name)

% This function performs a non-rigid FFD-based image registration
% to calculate the displacements within the timelapse

f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Displacement calculation';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');


%Access the files in the folder
folders = dir([path_names.store_path filesep experiment_name filesep folder_names.shift_corrected]);
folders = folders(3:end);

%Get the channels that need to be used
if not(strcmp(dispParam.marker_type,folder_names.fibers)) && ...
        not(strcmp(dispParam.marker_type,folder_names.beads))
    channel_names = {folder_names.beads,folder_names.fibers};
else
    channel_names = {dispParam.marker_type};
end

%Keep only the ones that are actually folders
tmp_indx_rmv = [];
for ii = 1:length(folders)
    %Remove the elements that are not folders, and that are not part of the
    %selected channels
    if not(folders(ii).isdir) || not(ismember(folders(ii).name,channel_names))
        tmp_indx_rmv = [tmp_indx_rmv ii];
    end
end
folders(tmp_indx_rmv) = [];
clear tmp_indx_rmv;

%Check if there is a cell mask to be used
if not(isfolder([path_names.store_path filesep experiment_name filesep folder_names.shift_corrected filesep folder_names.cellSeg]))
    dispParam.mask.cell = false;
end

% Parse the input/default options for the Parameter File
res = getResolution(folder_names.tiff_raw);
dispParam.imDim = length(res);
dispParam.resolution.xy = res(1);
if dispParam.imDim == 3
    dispParam.resolution.z = res(3);
end

dispParam.mask.flag = or(dispParam.mask.beads,dispParam.mask.cell);
dispParam.outImType = 'nii'; % alternatively 'tiff';

%Create the displacement results folder and store the parameter file
save_path = [path_names.store_path filesep experiment_name filesep folder_names.displacements];
if not(isfolder(save_path))
    mkdir(save_path);
end
save([save_path filesep disp_param_file_name],'dispParam');

%We save the grid size in microns, but we pass it to pixels for the
%operations
dispParam.grid.minSpacing.xy = round(dispParam.grid.minSpacing.xy / res(1)); %um to px
dispParam.grid.minSpacing.z = round(dispParam.grid.minSpacing.z / res(3)); %um to px
clear res;

%Prepare the parameters that will be used in the templates
[replaceStrIn,replaceStrOut,deleteLineKeyword] = parseElastixDispCalcParam(dispParam);

% Define the name of the parameter template
template_file_txt = [path_names.templateDispCalcFilePath filesep 'elastixDispCalcParam_template.txt'];

if dispParam.jacobMatrixFlag
    template_file_bat = [path_names.templateDispCalcFilePath filesep  'elastixDispCalcRun_withJacobian_template.bat'];
else
    template_file_bat = [path_names.templateDispCalcFilePath filesep  'elastixDispCalcRun_template.bat'];
end

%Definitions for the command line file to be run
fixed_im_file_name = 'fixedIm';
moving_im_file_name = 'movingIm';
fixed_mask_name = [];
moving_mask_name = [];
if dispParam.mask.flag
    fixed_mask_name = 'fixedMask';
    moving_mask_name = 'movingMask';
end
im_format = 'tiff'; % default format

% Define a base name for the result folder
res_folder_base_name = 'elastixOutput';

%Loop through the marker folders
for ii = 1:length(folders)
    
    %Get the folder name
    channel_name = folders(ii).name;
    
    %Check the time points in the folder
    files = dir([folders(ii).folder filesep channel_name]);
    files = files(3:end);
    
    %Compute the number of time points
    nTP = length(files);
    
    %Build the elastix results folders
    % Main folder
    calc_folder = [path_names.store_path filesep experiment_name filesep folder_names.elastix_calc filesep 'elastixDispCalc_' channel_name];
    if ~exist(calc_folder,'dir')
        mkdir(calc_folder)
    end
    common_run_file = [calc_folder filesep 'elastixDispCalcRun.bat'];
    
    % Folder for the images
    im_folder = [calc_folder filesep 'images'];
    if ~exist(im_folder,'dir')
        mkdir(im_folder)
    end
    
    % Folder for the parameters
    param_folder = [calc_folder filesep 'parameters'];
    if ~exist(param_folder,'dir')
        mkdir(param_folder)
    end
    
    %Create a .txt file modifying the template with the current settings
    out_file = [param_folder filesep 'elastixDispCalcParam.txt'];
    deleteTextLines(deleteLineKeyword,template_file_txt,out_file);
    replaceMultiString(replaceStrIn,replaceStrOut,out_file,out_file);
    
    %Create a .bat file modifying the template with the current settings
    out_file = [calc_folder filesep 'elastixDispCalcRun.bat'];
    [rmv_line,str_in,str_out] = setElastixBatFileStrings(path_names.elastixLibFolder,im_folder,param_folder,fixed_im_file_name,moving_im_file_name,im_format,fixed_mask_name,moving_mask_name,dispParam);
    deleteTextLines(rmv_line,template_file_bat,out_file);
    replaceMultiString(str_in,str_out,out_file,out_file);
    clear rmv_line str_in str_out out_file;   

    
    switch dispParam.strategy
        
        case 'direct'
            %Export the Fixed image
            exportElastixIm([folders(ii).folder filesep channel_name],'tp_R.mat',im_folder,fixed_im_file_name,im_format,dispParam,path_names,experiment_name,folder_names,fixed_mask_name);
            
            %Loop through the time points
            for tt = 1:nTP
                
                %Get the file name
                file_name = files(tt).name;
                
                %Make sure that we don't use the relaxed time point
                if not(strcmp(file_name,'tp_R.mat'))
                    
                    
                    %Update progress bar
                    d.Value = tt/nTP;
                    d.Message = sprintf(['Displacement calculation at ' file_name(1:end-4) '...']);    
                    
                    %Export the Moving image
                    exportElastixIm([folders(ii).folder filesep channel_name],file_name,im_folder,moving_im_file_name,im_format,dispParam,path_names,experiment_name,folder_names,moving_mask_name);
                    
                    % Generate the run (command line) file for the current timepoint
                    current_run_file = [calc_folder filesep 'elastixDispCalcRun_' file_name(1:end-4) '.bat'];
                    
                    %Create the result folder for this timepoint
                    current_res_folder = [calc_folder filesep res_folder_base_name '_' file_name(1:end-4)];
                    if exist(current_res_folder,'dir')
                        rmdir(current_res_folder,'s');
                    end
                    mkdir(current_res_folder);                    

                    rLine{1} = [];
                    strIn{1} = 'savePath';
                    strOut{1} = current_res_folder;
                    deleteTextLines(rLine,common_run_file,current_run_file);
                    clear rLine;
                    replaceMultiString(strIn,strOut,current_run_file,current_run_file);
                    clear strIn strOut;
                    
                    % Calculate the displacements
                    [dispField,registeredIm,jacobianMat] = ffdImReg(current_run_file,current_res_folder,dispParam);
                    
                    % Store the results
                    if not(isfolder([save_path filesep channel_name]))
                        mkdir([save_path filesep channel_name]);
                    end
                    save([save_path filesep channel_name filesep file_name],'dispField','-v7.3');
                    
                    if dispParam.jacobMatrixFlag
                        save([save_path filesep channel_name filesep file_name],'-append','registeredIm','jacobianMat','-v7.3');
                    else
                        save([save_path filesep channel_name filesep file_name],'-append','registeredIm','-v7.3');
                    end
                    clear dispField registeredIm jacobianMat current_res_folder;
                    
                end
            end 
            
            
        case 'sequential'
            %Loop through the time points
            for tt = 1:nTP
                
                %Get the file name
                file_name = files(tt).name;
                
                %Update progress bar
                d.Value = tt/nTP;
                d.Message = sprintf(['Displacement calculation at ' file_name(1:end-4) '...']);                
                
                %Export the Moving image (current time point)
                exportElastixIm([folders(ii).folder filesep channel_name],file_name,im_folder,moving_im_file_name,im_format,dispParam,path_names,experiment_name,folder_names,moving_mask_name);
                
                %Export the Fixed image (previous time point, considering the time step)
                if tt == 1
                    %If it is the first time point, we use itself as a fixed
                    %image
                    exportElastixIm([folders(ii).folder filesep channel_name],file_name,im_folder,fixed_im_file_name,im_format,dispParam,path_names,experiment_name,folder_names,fixed_mask_name);
                else
                    %If not, we use the previous one
                    exportElastixIm([folders(ii).folder filesep channel_name],files(tt - dispParam.seq.tStep).name,im_folder,fixed_im_file_name,im_format,dispParam,path_names,experiment_name,folder_names,fixed_mask_name);
                end
                
                % Generate the run (command line) file for the current timepoint
                current_run_file = [calc_folder filesep 'elastixDispCalcRun_' file_name(1:end-4) '.bat'];
                current_res_folder = [calc_folder filesep res_folder_base_name '_' file_name(1:end-4)];
                if exist(current_res_folder,'dir')
                    rmdir(current_res_folder,'s');
                end
                mkdir(current_res_folder);
                strIn{1} = 'savePath';
                strOut{1} = current_res_folder;
                replaceMultiString(strIn,strOut,common_run_file,current_run_file);
                clear strIn strOut;
                
                % Calculate the displacements
                [dispField,registeredIm,jacobianMat] = ffdImReg(current_run_file,current_res_folder,dispParam);
                
                % Store the results
                save([save_path filesep channel_name filesep file_name],'-append','dispField','registeredIm','jacobianMat','-v7.3');
                clear dispField registeredIm jacobianMat;
                
            end
            
    end
    
    % Delete the temporal folder containing the images
    try
        rmdir(im_folder,'s');
    catch
        disp('Impossible to delete the temporal folder used by elastix to store the images under analysis');
    end
    
end

% Close the dialog box
close(d);
delete(f);

end % end of the main function




function exportElastixIm(image_path,image_name,im_folder,im_file_name,im_format,dispParam,path_names,experiment_name,folder_names,mask_name)
% Load the fixed image
im = load([image_path filesep image_name]);
im = im.im_corr;

% Write the fixed image to a file
if exist([im_folder filesep im_file_name '.' im_format],'file')
    delete([im_folder filesep im_file_name '.' im_format])
end
writeIm2File(im,im_folder,im_file_name,im_format);
clear tmp;

% Calculate the mask and write it to a file (if required)
if dispParam.mask.flag
    mask = ones(size(im),'uint8');
    if dispParam.mask.beads
        mask = mask .* beadMaskCalc(im);
    end
    if dispParam.mask.cell
        %Load the mask image
        tmp = load([path_names.store_path filesep experiment_name filesep folder_names.shift_corrected filesep folder_names.cellSeg filesep image_name]);
        cell_mask = tmp.cell_seg;
        clear tmp;
        
        if not(isempty(cell_mask))
            mask = mask .* uint8((cell_mask==0));
        end
        clear cell_mask;
    end
    
    %Write it to a file
    if exist([im_folder filesep mask_name  '.'  im_format],'file')
        delete([im_folder filesep mask_name  '.'  im_format])
    end
    writeIm2File(mask,im_folder,mask_name,im_format);
    clear mask;
end
clear im;



end

