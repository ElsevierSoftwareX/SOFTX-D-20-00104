%    This file is part of cellMat.
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
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


function checkMetadata(rInfo,sInfo)

if (rInfo.xyResolution ~= sInfo.xyResolution)
    error('The xyResolution of bead images is different.');
end

if isfield(rInfo,'zResolution')&& isfield(sInfo,'zResolution')
    if (rInfo.zResolution ~= sInfo.zResolution)
        error('The zResolution of bead images is different.');
    end
end

if (rInfo.zNum>1)
    if xor(rInfo.reverseFlag,sInfo.reverseFlag)
        error('The reverse flag for the Zstack of bead images is different.');
    end
end

