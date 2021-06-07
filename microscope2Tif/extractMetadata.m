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

function info = extractMetadata(reader,seriesID)


metadata = reader.getMetadataStore(); 
% Select the 1st series
reader.setSeries(seriesID);
% Extract metadata 
info.imSize = [metadata.getPixelsSizeY(seriesID).getValue(), metadata.getPixelsSizeX(seriesID).getValue()];
info.nZplanes = metadata.getPixelsSizeZ(seriesID).getValue();
info.nCh = metadata.getPixelsSizeC(seriesID).getValue();
info.nTsteps =  metadata.getPixelsSizeT(seriesID).getValue();
% Spatial Resolution
pixelSizeX = metadata.getPixelsPhysicalSizeX(seriesID).value(ome.units.UNITS.MICROMETER).doubleValue(); % in um
pixelSizeY = metadata.getPixelsPhysicalSizeY(seriesID).value(ome.units.UNITS.MICROMETER).doubleValue(); % in um
if info.nZplanes>1
    pixelSizeZ = metadata.getPixelsPhysicalSizeZ(seriesID).value(ome.units.UNITS.MICROMETER).doubleValue(); % in um
    info.resolution = [pixelSizeX, pixelSizeY, pixelSizeZ];
    clear pixelSizeX pixelSizeY pixelSizeZ
else
    info.resolution = [pixelSizeX, pixelSizeY];
    clear pixelSizeX pixelSizeY
end
% Bit depth
info.bitDepth = reader.getBitsPerPixel();

clear metadata





end
