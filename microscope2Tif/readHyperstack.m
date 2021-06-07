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



function out = readHyperstack(reader,info,ch,zRange,tRange)

if length(tRange)==1
    tRange = [tRange,tRange];
end
if length(zRange)==1
    zRange = [zRange,zRange];
end

tVec = tRange(1):tRange(2);
zVec = zRange(1):zRange(2);

zNum = length(zVec);
tNum = length(tVec);

out = zeros([info.imSize,zNum,tNum],'uint8');


for tt=1:tNum
    for zz=1:zNum
        planeIndx = reader.getIndex(zVec(zz)-1,ch-1,tVec(tt)-1) + 1;
        tmp = double(bfGetPlane(reader,planeIndx));
        out(:,:,zz,tt) = uint8(round(255*(tmp/(2^info.bitDepth -1)))); % we convert the input images to 8 bits
        clear tmp planeIndx
    end
end

out = squeeze(out);

