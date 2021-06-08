
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

function imProc(timePoints,pathParam,shiftParam,filtParam,varargin)
% This function process the input images:
% 1) Filtering of Bead/Fiber/Cell images (if required)
% 2) Shift calculation (if required) 
% 3) Shift correction (if required)
% 4) Z-projection / Surface Detection (if required/needed)
% 5) Cell segmentation


% Chek if the cell images have to be processed
if nargin>4
    cellFlag = true;
    cellfiltParam = varargin{1};
    cellsegParam = varargin{2};
else
    cellFlag =false;
end


%% Filter the Images

if cellFlag
    timeSeriesImFilt(timePoints,pathParam,filtParam,cellfiltParam,cellsegParam);
else
    timeSeriesImFilt(timePoints,pathParam,filtParam);
end

%%  Shift Calculation

if strcmp(shiftParam.method,'pre') || strcmp(shiftParam.method,'both')
    shiftInfo = timeSeriesShiftCalc(timePoints,pathParam,shiftParam.pre);  
end

%% Shift Correction

if strcmp(shiftParam.method,'pre') || strcmp(shiftParam.method,'both')
    timeSeriesShiftCorr(timePoints,shiftInfo,cellFlag,pathParam)   
end

%% Dimension Reduction of the Bead/Fiber images 


%% Cell Segmentation

if cellFlag
    timeSeriesCellSeg(timePoints,pathParam,cellsegParam);
end
