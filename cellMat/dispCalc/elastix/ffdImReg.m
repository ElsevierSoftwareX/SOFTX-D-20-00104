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


function [dispField,registeredIm,jacobianMat] = ffdImReg(runFile,resFolder,dispParam,varargin)
% Launch elastix to perform the image registration and calculate the displacements

outImFormat = dispParam.outImType;
resolution = dispParam.resolution;
jacobMatrixFlag = dispParam.jacobMatrixFlag;

% Run the command line file that run elastix with the corresponding
% parameters and images
verboseFlag = false;
if verboseFlag % show in Matlab's command window the internal results of elastix
    status = system(['"' runFile  '"'],'-echo');
else
    [status,~] = system(['"' runFile '"']);
end
if (status ~= 0)
    error('Elastix-based calculation failed. Check provided parameters or input images.')
end


% Extract the scaling factor for the displacements used by elastix
if nargin>3
    transformParamFile = varargin{1};
else
    transformParamFile = 'TransformParameters.0.txt';
end
strIn{1} = '(Spacing';
paramLine = extractMultiLine(strIn,[resFolder filesep transformParamFile]); 
scale = regexp(paramLine{1},'\d+\.?\d*|-\d+\.?\d*|\.?\d*','match');  % Extract the numeric values from the previously read strings
scale = cellfun(@str2num, scale);

% Pass the results back to Matlab
jacobianMat = [];
switch outImFormat
    case 'nii'
        dispField = readDispNii([resFolder filesep 'dispField' filesep 'deformationField.nii'],scale,resolution);
        registeredIm = readImNii([resFolder filesep 'regIm' filesep 'result.nii']); 
        if jacobMatrixFlag
            jacobianMat = readJacNii([resFolder filesep 'jacobianMat' filesep 'fullSpatialJacobian.nii'],scale,resolution);
        end
    case 'tiff'
        dispField = readDispTiff([resFolder filesep 'dispField' filesep 'deformationField.tiff'],scale,resolution);
        registeredIm = read3Dim([resFolder filesep 'regIm' filesep 'result.tiff']); 
end
registeredIm = uint8(round(stretch(registeredIm)));


% Remove the results from the hard disk to free some space
try
    rmdir([resFolder filesep 'dispField'],'s')
    rmdir([resFolder filesep 'regIm'],'s')
    if jacobMatrixFlag
        rmdir([resFolder filesep 'jacobianMat'],'s')
    end
catch
    disp('Impossible to delete the temporal folder used by elastix to compute the displacements and the registered image')
end



