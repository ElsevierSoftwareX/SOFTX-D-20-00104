%    This file is part of TFMLAB.
%     Copyright (C) 2016-2020 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
%     Copyright (C) 2019-2020 Jorge Barrasa-Fano
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


function cellMask = cellSegmentation(cellIm,param)

%%%%%% param.thresholdMethod = 'otsu'; --> Thresholding using maximal inter-class variance method by Otsu (1979)
%%%%%% param.thresholdMethod = 'isodata'; --> Thresholding using the Isodata algorithm by Ridler and Calvard (1978)
%%%%%% param.thresholdTotNum --> A poisitve integer representing the total number of thresholds to compute if isodata algorithm is selected. Each threshold is identified by a positive integer between 1 and thresholdTotNum
%%%%%% param.thresholdLevel --> A positive integer (smaller than thresholdNum) representing the minimum threshold num to be used for thesegmentation of the cell. Only valid if (param.thresholdNum>1) and (param.thresholdMethod='isodata')
%%%%%% param.minCellSize --> A positive integer with the minimum number of pixels/voxels that a cell body should contain. If [] (empty), just the larger cell bpdy is kept.

try
    
    [~,thOtsu] = threshold(cellIm,'otsu');
    th = thOtsu - (thOtsu*(param.thPercentReduct/100));
    cellMask =threshold(cellIm,'fixed',th);    
    
    % Just keep the largest segmented object
    lb = label(cellMask);
    objSize = measure(lb,[],'size');
    if not(isempty(objSize))
        if isempty(param.minCellSize)
            [~,indx]=max(objSize.size);
            cellMask = lb==(objSize.id(indx));
        else
            indx=find(objSize.size>=param.minCellSize);
            cellMask = newim(size(lb),'bin');
            if not(isempty(indx))
                for ob=1:length(indx)
                    cellMask = cellMask + (lb==(objSize.id(indx(ob))));
                end
            end
        end
        clear lb indx
    end
    clear objSize
    % Fill holes within the mask
    cellMask = fillholes(cellMask);
    % Smooth out the mask shape
    if param.smoothCellMask
        gm = gradmag(cellMask,2);
        skl = bskeleton(threshold(gm));
        cellMask = fillholes(skl);
        clear gm skl
    end
    
catch
    
    s = size(cellIm);
    cellMask = newim(s,'bin') ; 
    
end



