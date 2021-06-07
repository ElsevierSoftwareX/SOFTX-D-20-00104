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

function im3D = read3Dim(fileName)
% This function reads multipage tiff images

warning('off','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning')
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')




mtdt = extractTifMetadata(fileName);


switch mtdt.bitDepth
    case 8
        dataType = 'uint8';
    case 12
        dataType = 'uint16';
    case 16
        dataType = 'uint16';
    case 24
        dataType = 'uint16';
end

im3D=zeros(mtdt.rowNum,mtdt.colNum,mtdt.imNum,dataType);

TifLink = Tiff(fileName, 'r');
for ii=1:mtdt.imNum
    TifLink.setDirectory(ii);
    im3D(:,:,ii)=TifLink.read();
end
TifLink.close();
clear  mtdt TifLink ii

warning('on','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning')
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')


