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



function info=fileName2metadata(fileName)

% Find out the dimension of the experiment
zFlag = false; tFlag = false;
indx = strfind(fileName,'XYZT');
if isempty(indx)
    indx = strfind(fileName,'XYTZ');
    if isempty(indx)
        indx = strfind(fileName,'XYZ');
        if isempty(indx)
            indx = strfind(fileName,'XYT');
            if isempty(indx)
                indx = strfind(fileName,'XY');
                hyperStack = 'XY';
            else
                tFlag = true;
                hyperStack = 'XYT';
            end
        else
            zFlag=true;
            hyperStack = 'XYZ';
        end
    else
        zFlag=true; tFlag =true;
        hyperStack = 'XYTZ';
    end   
else
    zFlag=true; tFlag =true;
    hyperStack = 'XYZT';
end
% Check if the Zstack has to be reversed
if zFlag
    indx = strfind(fileName,'rvs');
    if isempty(indx)
        info.reverseFlag = false;
    else
        info.reverseFlag = true;
    end
end
%
switch hyperStack
    case 'XYZT'
        info.hyperStack = 'ZT'; 
    case 'XYTZ'
        info.hyperStack = 'TZ'; 
    otherwise
        info.hyperStack = 'ZT';
end
info.zFlag = zFlag;
info.tFlag = tFlag;
clear indx hyperStack zFlag tFlag
% Extract resolution
xyresStartIndx = strfind(fileName,'xyres') + length('xyres');
info.xyResolution = str2num(fileName(xyresStartIndx:xyresStartIndx+3));
if info.zFlag
    zresStartIndx = strfind(fileName,'zres') + length('zres');  
    info.zResolution = str2num(fileName(zresStartIndx:zresStartIndx+3));
end
% Extract number of time points (if needed)
if (info.tFlag)&&(info.zFlag)
        tnumStartIndx = strfind(fileName,'TimePoints') + length('TimePoints');
        keywordVec = {'relaxed';'stressed';'cell';'relaxing'}; 
        for kk=1:length(keywordVec)
            if not(isempty(strfind(fileName(tnumStartIndx:end),keywordVec{kk})))
                keywordIndx = kk;
            end
        end
        tnumEndIndx = (tnumStartIndx-1) + strfind(fileName(tnumStartIndx:end),['_' keywordVec{keywordIndx}])-1;
        info.tNum = str2num(fileName(tnumStartIndx:tnumEndIndx));
        mtdt = extractTifMetadata(fileName);
        info.zNum = mtdt.imNum/info.tNum;
        clear tnumStartIndx tnumEndIndx mtdt
        
elseif (info.tFlag) && not(info.zFlag)
        mtdt = extractTifMetadata(fileName);
        info.tNum = mtdt.imNum;
        info.zNum = 1;
        clear mtdt
        
elseif not(info.tFlag) && (info.zFlag)
        mtdt = extractTifMetadata(fileName);
        info.zNum = mtdt.imNum;
        info.tNum = 1;
        clear mtdt
        
elseif not(info.tFlag) && not(info.zFlag)
        info.tNum =1;
        info.zNum =1;
end


