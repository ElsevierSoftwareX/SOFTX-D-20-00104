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

function mtdt = extractTifMetadata(fileName)


warning('off','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning')
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')


InfoImage=imfinfo(fileName);
mtdt.colNum=InfoImage(1).Width;
mtdt.rowNum=InfoImage(1).Height;
mtdt.imNum=length(InfoImage);
mtdt.bitDepth=InfoImage(1).BitDepth;


warning('on','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning')
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')