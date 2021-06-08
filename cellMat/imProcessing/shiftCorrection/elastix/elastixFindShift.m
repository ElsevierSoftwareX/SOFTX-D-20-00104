%    This file is part of cellMat.
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
%     Copyright (C) 2019 Jorge Barrasa Fano
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


function [shft,ax_ang] = elastixFindShift(imA,imB,options,file_name)

% This function writes the required parameters to run elastix and use it to
% estimate the existing translational shift betwen imA and imB



%%  Build required temporal folders for calculations

% Main temporal folder
tmpFolder = [options.calcFolder filesep 'elastixFindShift'];
if ~exist(tmpFolder,'dir')
    mkdir(tmpFolder)
end
% Folder for the images
imFolder = [tmpFolder filesep 'images'];
if ~exist(imFolder,'dir')
    mkdir(imFolder)
end
% Folder for the parameters
paramFolder = [tmpFolder filesep 'parameters'];
if ~exist(paramFolder,'dir')
    mkdir(paramFolder)
end
% Folder for the results
resFolder = [tmpFolder filesep 'elastixOutput_' file_name];
if exist(resFolder,'dir')
   rmdir(resFolder,'s')   
end
mkdir(resFolder)



%% Write the images to the corresponding folder

imFormat = 'tiff'; % default format
fixedImFileName = 'fixedIm';
movingImFileName = 'movingIm';
writeIm2File(imA,imFolder,fixedImFileName,imFormat);
writeIm2File(imB,imFolder,movingImFileName,imFormat);
fixedImFileName = [fixedImFileName '.' imFormat]; 
movingImFileName = [movingImFileName '.' imFormat];


%% Parameter file

% Parse the input/default options for the Parameter File
options.imDim = length(squeeze(size(imA)));
options.outImType = 'tiff';
[replaceStrIn,replaceStrOut,deleteLineKeyword] = parseElastixFindShiftParam(options);
% Copy and Modify the parameter file with the selected options
inFile = [options.templateFilePath filesep 'elastixFindShiftParam_template.txt'];
outFile = [paramFolder filesep 'elastixFindShiftParam.txt'];
deleteTextLines(deleteLineKeyword,inFile,outFile);
replaceMultiString(replaceStrIn,replaceStrOut,outFile,outFile);


%% Running (command line) file

% Copy and Modify the running (command line) file accordingly
inFile = [options.templateFilePath filesep 'elastixFindShiftRun_template.bat'];
outFile = [tmpFolder filesep 'elastixFindShiftRun.bat'];
strIn{1} = 'elastixLibPath'; strOut{1} = options.libFolder;
strIn{2} = 'imReadPath'; strOut{2} = imFolder;
strIn{3} = 'parametersReadPath'; strOut{3} = paramFolder;
strIn{4} = 'savePath'; strOut{4} = resFolder;
strIn{5} = 'fixedImName'; strOut{5} = fixedImFileName;
strIn{6} = 'movingImName'; strOut{6} = movingImFileName;
strIn{7} = 'parameterFileName'; strOut{7} = 'elastixFindShiftParam.txt';
replaceMultiString(strIn,strOut,inFile,outFile);
clear strIn strOut inFile outFile


%% Calculate the shift usng elastix
% Run the command line file that calls elastix with the required parameters

verboseFlag = false;
if verboseFlag % show in Matlab's command window the internal results of elastix
    status = system(['"' tmpFolder filesep 'elastixFindShiftRun.bat"'],'-echo');
else
    [status,~] = system(['"' tmpFolder filesep 'elastixFindShiftRun.bat"']);
end
if (status ~= 0)
    error('Elastix-based shift correction failed.Check provided parameters or input images.')
end


%% Import the results to Matlab

transformParamFile = 'TransformParameters.0.txt';
% Extract the line from the transformParamFile that contains the needed keywords
strIn{1} = '(TransformParameters';
strIn{2} = '(Spacing';
paramLine = extractMultiLine(strIn,[resFolder filesep transformParamFile]); 
% Extract the numeric values from the previously read strings
shft = regexp(paramLine{1},'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match');
scale = regexp(paramLine{2},'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match');  % scaling factor for the displacements used by elastix
shft = cellfun(@str2num, shft);
scale = cellfun(@str2num, scale);


%Rotations
if length(shft)>3
    rots = shft(1:3);
    shft = shft(4:end);
end

%Convert the displacements to pixels/voxels
shft = shft./scale;

%Make the displacements compatible with the sign convention used by the alternative methods
shft = -shft;

ax_ang = [];
if strcmp(options.txType,'EulerTransform')
    %Transform the Euler angles to rotation matrix
    rotm = eul2rotm(rots,"XYZ");    
    %Transform the rotation matrix into an axis angle
    ax_ang = rotm2axang(rotm);
end


clear strIn paramLine;













