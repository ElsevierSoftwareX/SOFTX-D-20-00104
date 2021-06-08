
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


function timeSeriesDispShiftCorr(path_names,experiment_name,folder_names,dispParam,file_names)
% This function corrects for any remaining rigid shift in the images

%The displacements, registeredIm and Jacobian matrices are overwritten in
%the displacements folder (where beads and fibers are separated, so there
%will be no conflicts between different shifts corrections for the
%different channels).

%The images tp_01.mat... are saved in a new post-shift correction folder
%according to their marker name.

%Create a waiting bar window
f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Displacement calculation (shift correction)';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');

%Access the files in the folder
folders = dir([path_names.store_path filesep experiment_name filesep folder_names.displacements]);
folders = folders(3:end);

%Create the displacement results folder and store the parameter file
save_path = [path_names.store_path filesep experiment_name filesep folder_names.post_shift_corrected];
if not(isfolder(save_path))
    mkdir(save_path);
end

%Keep only the ones that are actually folders
tmp_indx_rmv = [];
for ii = 1:length(folders)
    %Remove the elements that are not folders, and that are not part of the
    %selected channels
    if not(folders(ii).isdir)
        tmp_indx_rmv = [tmp_indx_rmv ii];
    end
end
folders(tmp_indx_rmv) = [];
clear tmp_indx_rmv;

%Loop through the channels (that are only the markers in fact)
for ii = 1:length(folders)
    
    %Get the folder name
    channel_name = folders(ii).name;
    
    %Check the time points in the folder
    files = dir([folders(ii).folder filesep channel_name]);
    files = files(3:end);
    
    %Compute the number of time points
    nTP = length(files);
    
    for tt = 1:nTP
        
        %Update progress bar
        d.Value = tt/nTP;
        d.Message = sprintf([channel_name ': calculating the shifts for time point ' num2str(tt) '...']);
        
        %Load the displacement field
        load([folders(ii).folder filesep channel_name filesep files(tt).name],'dispField');    
        
        %Compute the mean shift
        shft = [mean(dispField.X(:)),mean(dispField.Y(:))];
        if isfield(dispField,'Z')
            shft = [shft,mean(dispField.Z(:))] ;
        end
        clear dispField;
        
        %Initialize a container (if we are in the first timepoint). We do
        %it in this way since it's the most efficient way to already know
        %the dimension of the images
        if tt == 1
            shft_time_lapse = zeros(nTP,length(shft));  
        end
        shft_time_lapse(tt,:) = shft;
        clear shft;
    end   
    
    %Get resolution
    res = getResolution(folder_names.tiff_raw);
    
    %Convert to pixels
    shft_time_lapse_px = -shft_time_lapse ./ repmat(res,size(shft_time_lapse,1),1);
    
    %Save the shifts
    save([save_path filesep channel_name '_' file_names.shift_param],'shft_time_lapse_px');  
    
    %Compute the crops
    crop_start = zeros(size(shft_time_lapse_px));
    crop_end = zeros(size(shft_time_lapse_px));
    crop_start(shft_time_lapse_px>0) = ceil(shft_time_lapse_px(shft_time_lapse_px>0));
    crop_end(shft_time_lapse_px<0) = ceil(abs(shft_time_lapse_px(shft_time_lapse_px<0)));
    crop_startT = max(crop_start,[],1) ;
    crop_endT = max(crop_end,[],1);
    
    % The reference images does not have to be shifted
    ref_shift = zeros(1,size(shft_time_lapse_px,2));
    
    clear shft_time_lapse_px crop_start crop_end;
    
    %Loop through the displacement folder and shift correct it
    for tt = 1:nTP
        
        %Update progress bar
        d.Value = tt/nTP;
        d.Message = sprintf([channel_name ': correcting the shifts for the displacement field of time point ' num2str(tt) '...']);
        
        %Load the displacement field
        load([folders(ii).folder filesep channel_name filesep files(tt).name],'dispField','registeredIm','jacobianMat');
        
        switch dispParam.strategy
            case 'direct'
                %Correct the displacement field
                dispField.X = dispField.X - shft_time_lapse(tt,1);
                dispField.Y = dispField.Y - shft_time_lapse(tt,2);
                dispField.X = shiftCorr(dispField.X,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
                dispField.Y = shiftCorr(dispField.Y,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
                if isfield(dispField,'Z')
                    dispField.Z = dispField.Z - shft_time_lapse(tt,3);
                    dispField.Z = shiftCorr(dispField.Z,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
                end
                
                %Correct the registered Im
                registeredIm = shiftCorr(registeredIm,ref_shift,crop_startT,crop_endT,'im');
                
                %Check if the Jacobian Matrix exists
                if exist('jacobianMat','var') 
                    if not(isempty(jacobianMat))
                        %Correct the Jacobian Matrix
                        for jj = 1:numel(jacobianMat)
                            jacobianMat{jj} = shiftCorr(jacobianMat{jj},ref_shift,crop_startT,crop_endT,'field');
                        end
                    end
                        %Save everything back into the displacement folder
                        save([folders(ii).folder filesep channel_name filesep files(tt).name],'dispField','registeredIm','jacobianMat','-v7.3');                    
                else                    
                    %Save everything back into the displacement folder
                    save([folders(ii).folder filesep channel_name filesep files(tt).name],'dispField','registeredIm','-v7.3');
                end
            case 'sequential'
             %TO DO!!  
        end
        
        clear registeredIm jacobianMat dispField;
    end
end

%Access the files in the folder
folders = dir([path_names.store_path filesep experiment_name filesep folder_names.shift_corrected]);
folders = folders(3:end);

%Get the channels that need to be used
if not(strcmp(dispParam.marker_type,folder_names.fibers)) && ...
        not(strcmp(dispParam.marker_type,folder_names.beads))
    marker_names = {folder_names.beads,folder_names.fibers};
else
    marker_names = {dispParam.marker_type};
end

%Keep only the ones that are actually folders
tmp_indx_rmv = [];
for ii = 1:length(folders)
    %Remove the elements that are not folders, and that are not part of the
    %selected channels
    if not(folders(ii).isdir)
        tmp_indx_rmv = [tmp_indx_rmv ii];
    end
end
folders(tmp_indx_rmv) = [];
clear tmp_indx_rmv;

%Create the results folder for the channels
for ii = 1:length(marker_names)    
    tmp = [save_path filesep marker_names{ii}];
    if not(isfolder(tmp))
        mkdir(tmp);
    end
    clear tmp;
end



%Loop through the shift corrected channels
for ii = 1:length(folders)
    %Get the folder name
    channel_name = folders(ii).name;
    
    %Check the time points in the folder
    files = dir([folders(ii).folder filesep channel_name]);
    files = files(3:end);
    
    %Compute the number of time points
    nTP = length(files);
    
    %Loop through the time points
    for tt = 1:nTP
        %Make sure that the variables that we are loading do not exist in
        %the workspace (this might sound naive but since we are in a loop,
        %it could be important, since variables are called differently
        %depending on the channel (cell_seg / im_corr)).
        clear im_corr cell_seg;
        
        %Update progress bar
        d.Value = tt/nTP;
        d.Message = sprintf([channel_name ': correcting the shifts for the beads and cell images of time point ' num2str(tt) '...']);
        
        %Load the image 
        load([folders(ii).folder filesep channel_name filesep files(tt).name]);
        
        %Loop through the markers (since they have different shift
        %correction values)
        for mm = 1:length(marker_names)
            clear shft_time_lapse_px;
            
            %Get the marker name
            marker_name = marker_names{mm};
            
            %Load the shift information
            load([save_path filesep marker_name '_' file_names.shift_param]);
            
            %Compute the crops
            crop_start = zeros(size(shft_time_lapse_px));
            crop_end = zeros(size(shft_time_lapse_px));
            crop_start(shft_time_lapse_px>0) = ceil(shft_time_lapse_px(shft_time_lapse_px>0));
            crop_end(shft_time_lapse_px<0) = ceil(abs(shft_time_lapse_px(shft_time_lapse_px<0)));
            crop_startT = max(crop_start,[],1) ;
            crop_endT = max(crop_end,[],1);
            
            %Clear the remaining variables. We also clear im_corr and cell_seg
            %which are the names of the files within the shift corrected .mats.
            %We do this to make sure that the
            clear crop_start crop_end;
            
            %Create the output folder
            tmp = [save_path filesep marker_name filesep channel_name];
            if not(isfolder(tmp))
                mkdir(tmp);
            end
            
            %Check if its the relaxed time point
            if tt < nTP
                %Correct the image
                if exist('im_corr','var')
                    im_corr = shiftCorr(im_corr,shft_time_lapse_px(tt,:),crop_startT,crop_endT,'im');
                    save([tmp filesep files(tt).name] ,'im_corr');
                else
                    cell_seg = shiftCorr(cell_seg,shft_time_lapse_px(tt,:),crop_startT,crop_endT,'bin');
                    save([tmp filesep files(tt).name],'cell_seg');
                end
            else
                %This is the relaxed time point so no shifts are done, just
                %the crop
                
                %Correct the image
                if exist('im_corr','var')
                    im_corr = shiftCorr(im_corr,ref_shift,crop_startT,crop_endT,'im');
                    save([tmp filesep files(tt).name],'im_corr');
                else
                    cell_seg = shiftCorr(cell_seg,ref_shift,crop_startT,crop_endT,'bin');
                    save([tmp filesep files(tt).name],'cell_seg');
                end
            end            
        end
    end
end
% Close the dialog box
close(d);
delete(f);

%% DEPRECATED (2016)

% %% Shift Corr.
% 
% % We need to change the sign to be compatible with other methods
% 
% % load([pathParam.savePath filesep pathParam.saveSubFolder filesep pathParam.saveFileName '_tR.mat'],'dispParam')
% switch size(shft_time_lapse,2)
%     case 2
%         shftTimeLapsePx = -(1/dispParam.resolution.xy)*shft_time_lapse;
%     case 3
%         shftTimeLapsePx = -[(1/dispParam.resolution.xy)*shft_time_lapse(:,1), (1/dispParam.resolution.xy)*shft_time_lapse(:,2), (1/dispParam.resolution.z)*shft_time_lapse(:,3)];
% end
% clear dispParam    
% 
% % Compute the total required crop
% crop_start = zeros(size(shftTimeLapsePx));
% crop_end = zeros(size(shftTimeLapsePx));
% crop_start(shftTimeLapsePx>0) = ceil(shftTimeLapsePx(shftTimeLapsePx>0)); 
% crop_end(shftTimeLapsePx<0) = ceil(abs(shftTimeLapsePx(shftTimeLapsePx<0))); 
% crop_startT = max(crop_start,[],1) ;
% crop_endT = max(crop_end,[],1);
% 
% % Check if the shiftInfo variable exists within the .mat file
% w = whos('-file',[pathParam.savePath filesep pathParam.saveSubFolder filesep pathParam.saveFileName '_t' num2str(timePoints(1)) '.mat']);
% shftInfoFlag = false;
% for ii=1:length(w)
%     if strcmp(w(ii).name,'shiftInfo')
%         shftInfoFlag = true;
%     end
% end
% clear w
% 
% % The reference images do not have to be shifted
% ref_shift = zeros(1,size(shftTimeLapsePx,2));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STRESSED STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% for tt=1:tpNum
%     
%     tp = timePoints(tt);
%     currentFileName = [pathParam.saveFileName '_t' num2str(tp)];
%     
%     disp(' ... ... ...')
%     disp(['Disp Corr - Shift Corr : TIME POINT ' num2str(tp) ' ... ...'])
%     
%     
%     % Parse the shftTimeLapse variable to shiftInfo struct
%     if shftInfoFlag
%         load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'shiftInfo') 
%     end
%     shiftInfo.disp = shftTimeLapsePx;
% 
%       
%     % Correct the stressed image
%     load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'stressedIm')
%     stressedIm = shiftCorr(stressedIm,shiftInfo.disp(tt,:),crop_startT,crop_endT,'im'); % Correct the shift
%     save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','stressedIm','shiftInfo') % Store the results
%     clear stressedIm
%     
%     % Correct the registered image
%     switch dispCalcStrategy
%         case {'direct','bwdIniTf'}
%             load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'registeredIm')
%             registeredIm = shiftCorr(registeredIm,ref_shift,crop_startT,crop_endT,'im'); % Correct the shift
%             save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','registeredIm') % Store the results
%             clear registeredIm
%         case 'seqBackProp'
%             load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'registeredIm')
%             registeredIm = shiftCorr(registeredIm,ref_shift,crop_startT,crop_endT,'im'); % Correct the shift
%             save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','registeredIm') % Store the results
%             clear registeredIm
%             load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'registeredImSeq','dispParam')
%             fixtpindx = (tt+dispParam.seqBackProp.tStep); clear dispParam
%             if (fixtpindx)>tpNum
%                 registeredImSeq = shiftCorr(registeredImSeq,ref_shift,crop_startT,crop_endT,'im'); % Correct the shift
%             else
%                 registeredImSeq = shiftCorr(registeredImSeq,shiftInfo.disp(fixtpindx,:),crop_startT,crop_endT,'im'); % Correct the shift
%             end
%             save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','registeredImSeq') % Store the results
%             clear registeredImSeq fixtpindx
%     end
%     
%     % Cell image for display
%     load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'cellIm'); % Read the required image
%     if not(isempty(cellIm))
%         cellIm = shiftCorr(cellIm,shiftInfo.disp(tt,:),crop_startT,crop_endT,'im'); % Correct the shift
%         save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','cellIm') % Store the results
%     end
%     clear cellIm
%     % Cell image for segmentation
%     load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'cellSeg'); % Read the required image
%     if not(isempty(cellSeg))
%         cellSeg = shiftCorr(cellSeg,shiftInfo.disp(tt,:),crop_startT,crop_endT,'bin'); % Correct the shift
%         save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','cellSeg') % Store the results
%     end
%     clear cellSeg
%     
%     % Correct the displacements 
%     switch dispCalcStrategy
%         case {'direct','bwdIniTf'}
%              load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'dispField');
%              dispField.X = dispField.X - shft_time_lapse(tt,1); 
%              dispField.Y = dispField.Y - shft_time_lapse(tt,2); 
%              dispField.X = shiftCorr(dispField.X,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
%              dispField.Y = shiftCorr(dispField.Y,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
%              if isfield(dispField,'Z')
%                  dispField.Z = dispField.Z - shft_time_lapse(tt,3); 
%                  dispField.Z = shiftCorr(dispField.Z,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
%              end 
%              save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','dispField') % Store the results
%              clear dispField
%         case 'seqBackProp'
%              load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'dispField');
%              dispField.X = dispField.X - shft_time_lapse(tt,1); 
%              dispField.Y = dispField.Y - shft_time_lapse(tt,2); 
%              dispField.X = shiftCorr(dispField.X,-shiftInfo.disp(tt,:),crop_startT,crop_endT,'field'); % Correct the shift
%              dispField.Y = shiftCorr(dispField.Y,-shiftInfo.disp(tt,:),crop_startT,crop_endT,'field'); % Correct the shift
%              if isfield(dispField,'Z')
%                  dispField.Z = dispField.Z - shft_time_lapse(tt,3);
%                  dispField.Z = shiftCorr(dispField.Z,-shiftInfo.disp(tt,:),crop_startT,crop_endT,'field'); % Correct the shift
%              end 
%              save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','dispField') % Store the results
%              clear dispField
%              load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'dispFieldSeq');
%              dispFieldSeq.X = shiftCorr(dispFieldSeq.X,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
%              dispFieldSeq.Y = shiftCorr(dispFieldSeq.Y,ref_shift,crop_startT,crop_endT,'field'); % Correct the shift
%              if isfield(dispFieldSeq,'Z')
%                  dispFieldSeq.Z = shiftCorr(dispFieldSeq.Z,shiftInfo.disp(tt,:),crop_startT,crop_endT,'field'); % Correct the shift
%              end 
%              save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','dispFieldSeq') % Store the results
%              clear dispFieldSeq
%     end
%     
%     clear shiftInfo
%     
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RELAXED STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% refFile = [pathParam.savePath filesep pathParam.saveSubFolder filesep pathParam.saveFileName '_tR.mat'];
% if exist(refFile,'file')
%     
%     disp('Disp Corr - Shift Corrr : RELAXED ... ...')
%     currentFileName = [pathParam.saveFileName '_tR'];
%     
%     % Parse the shftTimeLapse variable to shiftInfo struct
%     if shftInfoFlag
%         load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'shiftInfo')
%     end
%     shiftInfo.disp = shft_time_lapse;
%     
%     
%     % Correct the relaxedIm too.
%     load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'relaxedIm'); % Read the required image
%     relaxedIm = shiftCorr(relaxedIm,ref_shift,crop_startT,crop_endT,'im'); % Correct the shift
%     save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','relaxedIm','shiftInfo') % Store the results
%     
%     % Cell image for display
%     load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'cellIm'); % Read the required image
%     if not(isempty(cellIm))
%         cellIm = shiftCorr(cellIm,ref_shift,crop_startT,crop_endT,'im'); % Correct the shift
%         save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','cellIm') % Store the results
%     end
%     clear cellIm
%     
%     % Cell image for segmentation
%     load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'cellSeg'); % Read the required image
%     if not(isempty(cellSeg))
%         cellSeg = shiftCorr(cellSeg,ref_shift,crop_startT,crop_endT,'bin'); % Correct the shift
%         save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','cellSeg') % Store the results
%     end
%     clear cellSeg
%     
% end
% 
% 
% %%%%%
% disp('ok')
