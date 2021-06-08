%    This file is part of TFMLAB.
%     Copyright (C) 2016-2020 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
%     Copyright (C) 2019-2020 Jorge Barrasa-Fano
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



function [shift_info, rot_info] = timeSeriesShiftCalc(path_names,experiment_name,folder_names,shiftParam,file_names)

% This function estimates the shift (rigid translation) within the timelapse

%If there is a fiber channel, it will use it for shift correction,
%otherwise it will use the beads.


f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Shift calculation';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');


if not(strcmp(shiftParam.method,'Translation phaseCorr (fast)')) 
    shiftParam.elastix.calcFolder = [path_names.store_path filesep experiment_name filesep folder_names.elastix_calc];
    shiftParam.elastix.templateFilePath = path_names.templateShiftCorrectionFilePath;
    shiftParam.elastix.libFolder = path_names.elastixLibFolder;
end

%Get the channel that we need to use for shift calculation
if(shiftParam.flag_fibers)
    channel_name = 'fibers';
else
    channel_name = 'beads';
end

%Access the files in the folder
files = dir([path_names.store_path filesep experiment_name filesep folder_names.tiff_filtered filesep channel_name]);
files = files(3:end);

%Choose a reference time point which is in the middle of the timelapse
nTP = length(files);
if nTP == 2
    ref_time_point_indx = 2;
else
    ref_time_point_indx = ceil(nTP/2);
end

%Create a folder to store the
%parameter file
save_path = [path_names.store_path filesep experiment_name filesep folder_names.shift_corrected];
if not(isfolder(save_path))
    mkdir(save_path);
end

switch shiftParam.strategy
    
    case  'seq' %Needs validation !! Needs implementation to handle rotations !!
        
        % Compute the sequential shifts
        for tt = 1:nTP
            
            %Update progress bar
            d.Value = tt/nTP;
            d.Message = sprintf(['Calculating shifts for ' files(tt).name(1:end-4) '...']);
            
            % Load the corresponding image A
            tt_previous = max(tt-1,1);
            imA = load([files(tt_previous).folder filesep files(tt_previous).name]);
            imA = imA.im;
            
            %Check if we are in the first timepoint
            if tt==1
                im_dim = length(size(imA));
                shift_info.seq = zeros(nTP,im_dim);
                rot_info.seq = []; % TO DO               
            else
                % Load the corresponding image B
                imB = load([files(tt).folder filesep files(tt).name]);
                imB = imB.im;
                % Compute the shift
                [shift_info.seq(tt,:),~] = shiftCalc(imA,imB,shiftParam,files(tt_previous).name(1:end-4));
            end
            clear imA imB;
        end
        % Compute the absolute shifts respect to the reference timepoint
        shift_info.ref = zeros(size(shift_info.seq));
        rot_info.ref = []; % TO DO
        if nTP == 2
            shift_info.ref(1,:) = - shift_info.seq(2,:);            
        else
            for tt = 1:(ref_time_point_indx-1)
                shift_info.ref(tt,:) = - sum(shift_info.seq(tt+1:ref_time_point_indx,:),1);
            end
            
            for tt=(ref_time_point_indx+1):nTP
                shift_info.ref(tt,:) = sum(shift_info.seq((ref_time_point_indx+1):tt,:),1);                
            end
        end
        
    case 'direct'
        %Load the reference
        imA = load([files(ref_time_point_indx).folder filesep files(ref_time_point_indx).name]);
        imA = imA.im;
        %Compute the dimension
        im_dim = length(size(imA));
        %Create a container for the shifts
        shift_info.ref = zeros(nTP,im_dim);
        rot_info.ref = zeros(nTP,im_dim+1);
        for tt = 1:nTP
             %Update progress bar
            d.Value = tt/nTP;
            d.Message = sprintf(['Calculating shifts for ' files(tt).name(1:end-4) '...']);
            
            if not(tt == ref_time_point_indx)
                % Load the corresponding image B
                imB = load([files(tt).folder filesep files(tt).name]);
                imB = imB.im;
                
                tic
                % Compute the shift
                [shft,ax_ang] = shiftCalc(imA,imB,shiftParam,files(tt).name(1:end-4));
                toc
                
                %Store the results
                shift_info.ref(tt,:) = shft;
                if not(isempty(ax_ang))
                    rot_info.ref(tt,:) = ax_ang;
                end
                %Save RAM
                clear  imB;
            end
        end
        %Save RAM
        clear imA;
end

%Include the reference timepoint indx
shift_info.ref_time_point_indx = ref_time_point_indx;

%Store the parameter file
save([save_path filesep file_names.shift_param],'shiftParam','shift_info');

% Close the dialog box
close(d);
delete(f);

end

