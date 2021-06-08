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

function dispField = readDispTiff(fileName,scale,res)


%% Read the metadata

fileInfo = bfGetReader(fileName); % to read metada
metadata = fileInfo.getMetadataStore(); % to read metada
fileInfo.setSeries(0); % Select the 0-th series
info.nIm =  fileInfo.getImageCount();
info.imSize = [metadata.getPixelsSizeY(0).getValue(), metadata.getPixelsSizeX(0).getValue()];
info.nZplanes = metadata.getPixelsSizeZ(0).getValue();
info.nTpoints = metadata.getPixelsSizeT(0).getValue(); % Sometimes Z-planes are wrongly interpreted as timepoints
info.nCh = metadata.getPixelsSizeC(0).getValue(); % this represent the dimension of the vecotr field (2 for X and Y components, and 3 if Z is also present);

%% Read the data


dispField.X = zeros([info.imSize,info.nZplanes,info.nTpoints]);
dispField.Y = zeros([info.imSize,info.nZplanes,info.nTpoints]);
if (info.nCh==3)
    dispField.Z = zeros([info.imSize,info.nZplanes,info.nTpoints]);
end
for zz=1:info.nZplanes
    for tt=1:info.nTpoints
        dispField.X(:,:,zz,tt) = readImages(fileInfo,zz,tt,1);
        dispField.Y(:,:,zz,tt) = readImages(fileInfo,zz,tt,2);
        if (info.nCh==3)
            dispField.Z(:,:,zz,tt) = readImages(fileInfo,zz,tt,3);
        end
    end
end
dispField.X = squeeze(dispField.X);
dispField.Y = squeeze(dispField.Y);
if (info.nCh==3)
    dispField.Z = squeeze(dispField.Z);
end

% close the reader
fileInfo.close();
clear fileInfo metadata

%% Convert the displacements to the right units

% Convert the displacements to pixels/voxels
dispField.X = dispField.X / scale(1);
dispField.Y = dispField.Y / scale(2);
if (info.nCh==3)
    dispField.Z = dispField.Z / scale(3);
end
% Finally convert them to microns
dispField.X = dispField.X*res.xy;
dispField.Y = dispField.Y*res.xy;
if (info.nCh==3)
    dispField.Z = dispField.Z*res.z;
end

end




function im = readImages(fileInfo,zPlane,timePoint,ch)

imIndx = loci.formats.FormatTools.getIndex(fileInfo, zPlane-1, ch-1, timePoint-1) + 1;
im = double(bfGetPlane(fileInfo, imIndx));

end
