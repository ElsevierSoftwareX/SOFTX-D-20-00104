
%    This file is part of cellMat.
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
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

function [gParam,flag,refTimePoint] =  imfile2mat(flag,pathParam,gParam,shiftParam,filtParam,cellfiltParam,cellsegParam,refTimePoint)
% This function reads the images for the whole timelapse and store them
% into .mat files
% If the metadata is extracted from the file name of the images, the values
% provided by the user under the struct variable gParam are updated.


%% Extract info

% Check if the Relaxed image has been provided
 if not(isempty(pathParam.readPath.relaxed)) && not(isempty(pathParam.readFile.relaxed))
     flag.relaxedIm = true;
 else
     flag.relaxedIm = false;
 end
% Check if the Cell Image has been provided
 if not(isempty(pathParam.readPath.cell)) && not(isempty(pathParam.readFile.cell))
     flag.cellIm = true;
 else
     flag.cellIm = false;
 end
 
if gParam.readFileNameInfo % Extract some metadata from file names

    sInfo =  fileName2metadata([pathParam.readPath.stressed filesep pathParam.readFile.stressed ]); % Stressed Beads
    if flag.relaxedIm
        rInfo = fileName2metadata([pathParam.readPath.relaxed filesep pathParam.readFile.relaxed ]); % Relaxed Beads
        checkMetadata(rInfo,sInfo);    % Check if fInfo matches sInfo
    end
    gParam.xyResolution = sInfo.xyResolution;
    if sInfo.zNum>1
        gParam.zResolution = sInfo.zResolution;
    else
        gParam.zResolution = [];
    end
    gParam.timepointNum = sInfo.tNum;
    if isfield(sInfo,'reverseFlag')
        gParam.reverseZstack = sInfo.reverseFlag;
    else
        gParam.reverseZstack = [];
    end
    if isfield(sInfo,'hyperStack')
        gParam.hyperstackOrder = sInfo.hyperStack;
    else
        gParam.hyperstackOrder = [];
    end
    
else % Get the information directly from the user provided parameters
    
    % Stressed Image
    sInfo.tNum = gParam.timepointNum;
    mtdt = extractTifMetadata([pathParam.readPath.stressed filesep pathParam.readFile.stressed ]);
    sInfo.zNum = mtdt.imNum/sInfo.tNum; clear mtdt
    if (sInfo.tNum>1)&&(sInfo.zNum>1)
        sInfo.hyperStack = gParam.hyperstackOrder;
    end
    sInfo.reverseFlag = gParam.reverseZstack;
    % Relaxed image
    if flag.relaxedIm
        rInfo.tNum = 1;
        mtdt = extractTifMetadata([pathParam.readPath.relaxed filesep pathParam.readFile.relaxed ]);
        rInfo.zNum = mtdt.imNum; clear mtdt
        rInfo.reverseFlag = gParam.reverseZstack;
    end
 
end


 % Check the timepoints that will be analyzed
if isempty(gParam.timepoints)
    gParam.timepoints = 1:sInfo.tNum;
end


%% Read the stressed images (and cell images if provided) and save them into matlab .mat files


readInfo.zNum = sInfo.zNum;
readInfo.tNum = sInfo.tNum;
readInfo.z = []; % empty implies reading all
if isfield(sInfo,'hyperStack')
    readInfo.hyperStack = sInfo.hyperStack;
end

disp(' ... ... ...')
disp(' ... ... ...')
    
for tt=1:length(gParam.timepoints)
    
    tp = gParam.timepoints(tt);
    saveFileName = [pathParam.saveFileName '_t' num2str(tp)];

    disp(' ... ... ...')
    disp(['Importing to Matlab: TIME POINT ' num2str(tp) '/' num2str(sInfo.tNum) ' ... ...'])
    
    % Read the stressed image for the current time point
    disp(['Reading the Stressed image for T=' num2str(tp) '...'])
    readInfo.t = tp;
    stressedIm = squeeze(readNDim([pathParam.readPath.stressed filesep pathParam.readFile.stressed],readInfo));
    % Reverse the Z-stack if required
    if (length(size(stressedIm))==3) && sInfo.reverseFlag
        stressedIm = flipdim(stressedIm,3);
    end
       
    % Read the cell image for the current time point
    if flag.cellIm 
        disp(['Reading the Cell image for T=' num2str(tp) '...'])
        readInfo.t = tp;
        cellIm = squeeze(readNDim([pathParam.readPath.cell filesep pathParam.readFile.cell],readInfo));
        % Reverse the Z-stack if required
        if (length(size(cellIm))==3) && sInfo.reverseFlag
            cellIm = flipdim(cellIm,3);
        end
    else
        cellIm = [];
    end
    cellSeg = [];
    
    % Store the the results
    save([pathParam.savePath filesep pathParam.saveSubFolder filesep saveFileName '.mat'],'-v7.3','flag','cellIm','cellSeg','stressedIm','pathParam','gParam','shiftParam','filtParam','cellfiltParam','cellsegParam')
    clear cellIm stressedIm
    
end


%% Read the relaxed image and save them into matlab .mat file

if isempty(refTimePoint)
    refCalc_flag = true;
else
    if isnan(refTimePoint)
        refCalc_flag = false; % if refTimePoint=NaN --> sequential displacements will be computed and an absolute reference (relaxed state) is not needed
    else
        refCalc_flag = true;
    end
end

if refCalc_flag % if refTimePoint=NaN --> sequential displacements will be computed and an absolute reference (relaxed state) is not needed
    
    saveFileName = [pathParam.saveFileName '_tR'];
    
    if flag.relaxedIm
        
        readInfo.zNum = rInfo.zNum;
        readInfo.tNum = rInfo.tNum;
        readInfo.z = []; % empty implies reading all
        
        % Read the relaxed image for the current time point
        disp('Reading the Relaxed image ...')
        readInfo.t = 1;
        relaxedIm = squeeze(readNDim([pathParam.readPath.relaxed filesep pathParam.readFile.relaxed],readInfo));
        % Reverse the Z-stack if required
        if (length(size(relaxedIm))==3) && rInfo.reverseFlag
            relaxedIm = flipdim(relaxedIm,3);
        end
        % No Cell Image
        cellIm = [];
        cellSeg = [];
        
        
    else % we use the timepoint specified by the user as the reference (the relaxed state)
        
        readInfo.zNum = sInfo.zNum;
        readInfo.tNum = sInfo.tNum;
        readInfo.z = []; % empty implies reading all
        if isfield(sInfo,'hyperStack')
            readInfo.hyperStack = sInfo.hyperStack;
        end
        
        % Read the timepoint specified by the user
        disp('Reading the Relaxed image ...')
        if isempty(refTimePoint)
            readInfo.t = sInfo.tNum; % last timepoint will be used as the relaxed state
            refTimePoint = sInfo.tNum;
        else
            readInfo.t = refTimePoint;
        end
        relaxedIm = squeeze(readNDim([pathParam.readPath.stressed filesep pathParam.readFile.stressed],readInfo));
        % Reverse the Z-stack if required
        if (length(size(relaxedIm))==3) && sInfo.reverseFlag
            relaxedIm = flipdim(relaxedIm,3);
        end
        
        % Read the cell image from the timepoint specified by the user (the reference)
        if flag.cellIm
            disp('Reading the Cell image for the Realxed state')
            if isempty(refTimePoint)
                readInfo.t = sInfo.tNum; % last timepoint will be used
            else
                readInfo.t = refTimePoint;
            end
            cellIm = squeeze(readNDim([pathParam.readPath.cell filesep pathParam.readFile.cell],readInfo));
            % Reverse the Z-stack if required
            if (length(size(cellIm))==3) && sInfo.reverseFlag
                cellIm = flipdim(cellIm,3);
            end
        else
            cellIm = [];
        end
        cellSeg = [];
        
    end
    
    % Store the the results
    save([pathParam.savePath filesep pathParam.saveSubFolder filesep saveFileName '.mat'],'-v7.3','flag','cellIm','cellSeg','relaxedIm','pathParam','gParam','shiftParam','filtParam','cellfiltParam','cellsegParam')
    clear relaxedIm
    
end % end of --> if refTimePoint~=NaN


%%

disp('ok')
