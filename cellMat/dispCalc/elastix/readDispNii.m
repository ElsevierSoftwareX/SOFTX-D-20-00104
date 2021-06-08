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

function dispField = readDispNii(fileName,scale,res)


%% Read nii data 

niiData = load_nii(fileName);
data = squeeze(niiData.img); clear niiData

dim = length(size(data));

%% Extract the displacements

switch dim
    case 3 % 2D
        dispField.X = squeeze(data(:,:,1));
        dispField.Y = squeeze(data(:,:,2));
    case 4 % 3D
        dispField.X = squeeze(data(:,:,:,1));
        dispField.Y = squeeze(data(:,:,:,2));
        dispField.Z = squeeze(data(:,:,:,3));
end

% Rotate approprietely the XY-plane 
dispField.X = rot90(dispField.X,-1);
dispField.Y = rot90(dispField.Y,-1);
if dim==4
    dispField.Z = rot90(dispField.Z,-1);
end

% Flip the rows
dispField.X = flip(dispField.X,1);
dispField.Y = flip(dispField.Y,1);
if dim==4
    dispField.Z = flip(dispField.Z,1);
end


%% Convert the displacements to the right units


% Convert the displacements to pixels/voxels
dispField.X = dispField.X / scale(1);
dispField.Y = dispField.Y / scale(2);
if dim==4
    dispField.Z = dispField.Z / scale(3);
end

% Finally convert them to microns
dispField.X = dispField.X*res.xy;
dispField.Y = dispField.Y*res.xy;
if dim==4
    dispField.Z = dispField.Z*res.z;
end




