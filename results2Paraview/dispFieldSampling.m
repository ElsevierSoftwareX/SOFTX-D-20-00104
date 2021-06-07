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

function indx = dispFieldSampling(fieldSize,param)


switch param.type
    case 'grid'
        indx = false(fieldSize);
        indx(1:param.grid.xy:end,1:param.grid.xy:end,1:param.grid.z:end) = true;
    case 'random'
        pointNum = round((param.random.percent/100)*prod(fieldSize));
        randIndx = randperm(prod(fieldSize));
        indx = randIndx(1:min(pointNum,prod(fieldSize)));
        clear randIndx pointNum
    case 'beads'
        tmp = load(param.beads.file,param.beads.varName);
        beadsIm = cropData(tmp.(param.beads.varName),param.beads.roi); clear tmp
        beadsIm = padarray(beadsIm,param.beads.padSize*ones(1,3),'replicate'); 
        beadPos = imregionalmax(beadsIm); clear beadsIm
        info = regionprops('table',beadPos,'centroid'); clear beadPos
        subIndx = round(info.Centroid); % [x,y,z]
        indx = sub2ind(fieldSize,subIndx(:,2),subIndx(:,1),subIndx(:,3)); clear info subIndx
end