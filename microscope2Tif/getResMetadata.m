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



function metadataStr = getResMetadata(info)


% Convert the xy resolution to a string
strpxXY = num2str(round(100*info.resolution(1))/100); % round the resolution to 2 decimals
if length(strpxXY)<4
    if length(strpxXY)==1
        strpxXY = [strpxXY '.' repmat('0',1,2)];
    else
        strpxXY = [strpxXY repmat('0',1,4-length(strpxXY))];
    end
elseif length(strpxXY)>4
    strpxXY = strpxXY(1:4);
end
% Convert the z resolution to a string
if info.nZplanes>1
    strpxZ = num2str(round(100*info.resolution(3))/100); % round the resolution to 2 decimals
    if length(strpxZ)<4
        if length(strpxZ)==1
            strpxZ = [strpxZ '.' repmat('0',1,2)];
        else
            strpxZ = [strpxZ repmat('0',1,4-length(strpxZ))];
        end
    elseif length(strpxZ)>4
        strpxZ = strpxZ(1:4);
    end
end

% Write useful metadata to string
if info.nZplanes>1
    metadataStr =  ['_xyres' strpxXY '_zres' strpxZ];
else
    metadataStr =  ['_xyres' strpxXY];
end


end