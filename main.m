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

clear;
close all;
clc;

%Run the deafult parameters
parameterFile;

%Add the guis to the path
addpath(genpath([pwd filesep code_names.guis]));
addpath(genpath([pwd filesep code_names.images]));
addpath(genpath([pwd filesep code_names.external]));
addpath(genpath([pwd filesep code_names.utilities]));


%Add other external software 
% addpath(genpath('C:\Workdir\Programs'));

%% FILE SELECTION (GUI)

%Create a flag to export to multitiff or not
flag_multitiff = false;

%Ask for file selection
h_gui = set_inputs(folder_names);
h_gui.UIFigure.Name = 'TFMLAB: Input data';
disp('Data selection: waiting for user inputs...')
waitfor(h_gui);

%Check if we have to close the app
if flag_closeApp
    return;
end
%% Conversion to Multitiff

if flag_multitiff
    addpath(genpath([pwd filesep code_names.microscope2Tif]));
    %Convert from .lif format to .tif
    folder_names = microscope2Tif(file_path,file_name,store_path,folder_names,param,file_names);
    rmpath(genpath([pwd filesep code_names.microscope2Tif]));
end
%% IMAGE PROCESSING PARAMETER SELECTION (GUI)
addpath(genpath([pwd filesep code_names.cellMat]));

%Filtering paramter tunning
if exist('store_path','var')
    %The user is following the regular workflow. They completed the
    %miscroscope2Tif section and they want to keep going forward with the
    %initial files
    
    %Start the selection of filter parameters
    h_gui = set_imProcessing(store_path,folder_names,param.experiment_name);
else
    %The user skipped the miscroscope2Tif section and wants to select
    %existing raw tiff files
    
    %Ask the user to select a folder where the raw tiffs are stored
    tmp = uigetdir('', ['Please select your ', folder_names.tiff_raw, ' folder']);
    
    %Check if the user closed the app
    if tmp == 0
        return;
    end
    
    %Get the experiment name and the store path
    [param.experiment_name,path_names.store_path,folder_names.tiff_raw] = getInfoFromSelectedFolder(tmp);
    clear tmp;
    
    %Start the selection of filter parameters
    h_gui = set_imProcessing(path_names.store_path,folder_names,param.experiment_name);
end

disp('Image filtering and segmentation: waiting for user inputs...')
h_gui.UIFigure.Name = 'TFM LAB: Filter parameters';
waitfor(h_gui);

if exist('store_path','var')
    %Save the store path in the official container
    path_names.store_path = store_path;
    clear store_path;
end

%Check if the user closed the app
if exist('need_to_load','var')
    %Now three situations may have happened:
    %1 - The user skipped the step because they already have some filtered images
    %2 - The user wants to use predefined filter parameters and then start filtering
    %3 - The user selected filtering parameters using the GUI and wants to start filtering
    
    switch need_to_load
        case 'images'
            %Situation 1
            
            %Ask the user for the folder where their filtered images are stored
            tmp = uigetdir([path_names.store_path filesep param.experiment_name], ['Please select your ', folder_names.tiff_filtered, ' folder']);
            %Check if the user closed the app
            if tmp == 0
                return;
            end
            %Get all the required fields
            [param.experiment_name,path_names.store_path,folder_names.tiff_filtered] = getInfoFromSelectedFolder(tmp);
            clear tmp;
            
        case 'parameters'
            %Situation 2
            
            %Ask the user for the file where their filter parameters are stored
            [tmp_file,tmp_path] = uigetfile([path_names.store_path filesep param.experiment_name], ['Please select your ', file_names.filt_param, ' file (normally at the ', folder_names.tiff_filtered, ' folder)']);
            
            %Load the filter paramters
            load([tmp_path filesep tmp_file]);
            clear tmp_file tmp_path;
    end
else
    return;
end

%Perform the filtering
if not(strcmp(need_to_load,'images'))
    tic
    timeSeriesImFilt(path_names.store_path,param.experiment_name,folder_names,file_names,filtParam,cellsegParam,cellfiltParam);
    disp(['Time for image filtering = ', num2str(toc)]);
end
clear need_to_load;
%% Shift correction

%Define a flag to determine when the user is satisfied with the shift
%correction
flag_need_shift_correction = true;

while flag_need_shift_correction
    %Start GUI for shift correction (to get the strategy and the method)
    h_gui = set_shiftCorrection(path_names.store_path,folder_names,param.experiment_name);
    h_gui.UIFigure.Name = 'TFMLAB: Shift correction parameters';
    waitfor(h_gui);
    
    if exist('strategy','var')
        
        % The user didn't skip the step
        shiftParam.strategy = strategy; % 'seq' or 'direct'
        shiftParam.method = method; %
        shiftParam.flag_fibers = flag_fibers;
        
        tic
        % Compute the shifts for the timelapse
        [shift_info, rot_info] = timeSeriesShiftCalc(path_names,param.experiment_name,folder_names,shiftParam,file_names);
        disp(['Time for shift calculation ', num2str(toc)]);
        
        tic
        % Correct the shifts for the timelapse
        timeSeriesShiftCorr(path_names,param.experiment_name,folder_names,shift_info,rot_info);
        disp(['Time for shift correction ', num2str(toc)]);
        
        
        %Start GUI for shift correction check
        h_gui = check_shiftCorrection(path_names.store_path,folder_names,param.experiment_name);
        h_gui.UIFigure.Name = 'TFM LAB: Shift correction check';
        waitfor(h_gui);
        
        
    else
        %The user skipped the step because they already have some shift
        %corrected images
        
        %Ask the user for the folder where their filtered images are stored
        tmp = uigetdir([path_names.store_path filesep param.experiment_name], ['Please select your ', folder_names.shift_corrected, ' folder']);
        %Check if the user closed the app
        if tmp == 0
            return;
        end
        
        %Get all the required fields
        [param.experiment_name,path_names.store_path,folder_names.shift_corrected] = getInfoFromSelectedFolder(tmp);
        clear tmp;
        flag_need_shift_correction = false;
    end
    
    if isnan(flag_need_shift_correction)
        %The user cancelled the program
        return;
    end
end
%% Cell segmentation
%Check if cell segmentation is needed
if not(isfolder([path_names.store_path filesep param.experiment_name filesep folder_names.shift_corrected filesep 'cellSeg'])) && ...
        isfolder([path_names.store_path filesep param.experiment_name filesep folder_names.shift_corrected filesep 'cell'])
    
    %Cell Segmentation is needed because there is a 'cell' channel but not
    %a 'cellSeg' channel
    
    if not(exist('cellsegParam','var'))
        %Load the cell segmentation parameters
        load([path_names.store_path filesep param.experiment_name filesep folder_names.tiff_filtered filesep file_names.filt_param],'cellsegParam');
    end
    tic
    %Perform the cell segmentation
    timeSeriesCellSeg(path_names,param.experiment_name,folder_names,cellsegParam);
    disp(['Cell segmentation ',num2str(toc)]);
end
%% Displacement calculation

%Start the selection of displacement calculation parameters
h_gui = set_dispCalc(path_names.store_path, folder_names, param.experiment_name,dispParam);

disp('Displacement calculation: waiting for user inputs...')
h_gui.UIFigure.Name = 'TFM LAB: Displacement calculation parameters';
waitfor(h_gui);

%Now three situations may have happened:
%1 - The user skipped the step because they already have some computed
%displacements
%2 - The user wants to use predefined displacement parameters and then start computing displacements
%3 - The user selected parameters using the GUI and wants to start computing displacements

switch need_to_load
    case 'images'
        %Situation 1
        
        %Ask the user for the folder where their filtered images are stored
        tmp = uigetdir([path_names.store_path filesep param.experiment_name], ['Please select your ', folder_names.displacements, ' folder']);
        %Check if the user closed the app
        if tmp == 0
            return;
        end
        %Get all the required fields
        [param.experiment_name,path_names.store_path,folder_names.displacements] = getInfoFromSelectedFolder(tmp);
        clear tmp;
        
    case 'parameters'
        %Situation 2
        
        %Ask the user for the file where their filter parameters are stored
        [tmp_file,tmp_path] = uigetfile([path_names.store_path filesep param.experiment_name], ['Please select your ', file_names.disp_param, ' file (normally at the ', folder_names.displacements, ' folder)']);
        
        %Load the filter paramters
        load([tmp_path filesep tmp_file]);
        clear tmp_file tmp_path;
end

%Perform the displacement calculation
if not(strcmp(need_to_load,'images'))
    tic
    timeSeriesDispCalc(path_names,param.experiment_name,folder_names,dispParam,file_names.disp_param);
    disp(['Displacement measurement ', num2str(toc)]);
    
    tic
    timeSeriesDispShiftCorr(path_names,param.experiment_name,folder_names,dispParam,file_names);
    disp(['Displacement shift correction: ' num2str(toc)]);
    
end
clear need_to_load;

rmpath(genpath([pwd filesep code_names.cellMat]));
%% Export results to Paraview
%Ask if the user wants to use existing crop
question_initial = 'Do you want to use an existing crop and go directly to traction calculation?';
question_title = 'Crop selection';
msg_negative = 'No, I want to create a new crop';
answer = questdlg(compose(question_initial),question_title,'Yes',msg_negative,'Cancel','Cancel');

 addpath(genpath([pwd filesep code_names.results2Paraview]));
switch answer
    
    case 'Yes'  %The user wants to use existing files
        
        [~,tmp]= uigetfile([path_names.store_path filesep param.experiment_name], ['Please select your ', file_names.paraview_param, ' file.'...
            '( It should be in the subfolders within ', folder_names.results_for_paraview, ' )']);
        
        %Check if the user closed the app
        if tmp == 0
            return;
        else
            %Save it in the paraview results folder (the same file for all
            %the channels
            channels = dir([path_names.store_path filesep param.experiment_name filesep folder_names.results_for_paraview]);
            channels = channels(3:end);
            for ii = 1:length(channels)
                if channels(ii).isdir
                    if not(exist([channels(ii).folder filesep channels(ii).name filesep file_names.paraview_param],'file'))
                        copyfile([tmp filesep file_names.paraview_param],[channels(ii).folder filesep channels(ii).name]);
                    end
                end
            end
        end        
        clear tmp channels ii;
        
    case msg_negative  
        
        results2paraview(path_names,param.experiment_name,folder_names,file_names,code_names);
        
        %Segment nuclei if there is any
%         if isfolder([path_names.store_path filesep param.experiment_name filesep folder_names.shift_corrected filesep folder_names.nuclei])
%             
%             %Open GUI for nuclei segmentation
%             h_gui = set_nucleiSegmentation(path_names.store_path, param.experiment_name, folder_names, file_names);
%             disp('Nuclei segmentation: waiting for user inputs...')
%             h_gui.UIFigure.Name = 'TFM LAB: Nuclei segmentation parameters';
%             waitfor(h_gui);
%         end
        
    otherwise
        return;
end

 rmpath(genpath([pwd filesep code_names.results2Paraview]));
%% Traction recovery

%First check for the visual studio path
[abq1_text,abq2_text] = checkCompatibility();

if not(isempty(abq1_text)) && not(isempty(abq2_text))
    addpath(genpath([pwd filesep code_names.tractionRecovery]));
    
    %First check if there are any cells
    if isfolder([path_names.store_path filesep param.experiment_name filesep folder_names.shift_corrected filesep folder_names.cell])
        
%         %Start the selection of traction calculation parameters
%         load([path_names.store_path filesep param.experiment_name filesep folder_names.results_for_paraview filesep file_names.paraview_param],'crop_rect');
%         
%         %Check how many channels are there
%         channels = dir([path_names.store_path filesep param.experiment_name filesep folder_names.post_shift_corrected]);
%         channels = channels(3:end);
%         tmp_ids = [channels(:).isdir];
%         channels = channels(tmp_ids);
%         clear tmp_ids cell_seg;
%         
%         %Create a variable to store the cell segmentations of each channel
%         cell_seg.beads = [];
%         cell_seg.fibers = [];
%         for ii=1:length(channels)
%             tmp = load([path_names.store_path filesep param.experiment_name filesep folder_names.post_shift_corrected filesep channels(ii).name filesep folder_names.cellSeg filesep 'tp_01.mat'],'cell_seg');
%             %Crop and store the cell segmentation
%             cell_seg.(channels(ii).name) = tmp.cell_seg(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
%             clear tmp;
%         end
        
        %Open the GUI
        h_gui = set_tracCalc(path_names,folder_names,param,file_names);
        clear channels cell_seg;
        
        disp('Traction calculation: waiting for user inputs...')
        h_gui.UIFigure.Name = 'TFM LAB: Traction calculation parameters';
        waitfor(h_gui);
        
%         mesh_options.quality = 1.414;
%         mesh_options.maxsize_cell=2;
%         mesh_options.maxvol_gel=200;
        tstart = tic;
        tractionRecovery(path_names.store_path, folder_names, file_names,param.experiment_name,mech_props,mesh_options,abq1_text,abq2_text,code_names);
        disp(['Force calculation: ' num2str(toc(tstart))]);
        
    end
    rmpath(genpath([pwd filesep code_names.tractionRecovery]));
end
%% Remove the functions from the path
rmpath(genpath(pwd));