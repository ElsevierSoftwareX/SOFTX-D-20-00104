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


function[filtCell,smoothVal] = cellFilter(cellIm,lightenDarkRegions,smoothVal,lb,ub)
% This function applies filtering, contrast stretching and Z-projection to the acquired cell image

% This is the workflow followed:
% 1- Log/root Mapping to lighten up dark regions (if required)
% 2- Smooth out the cell image (if required)
% 3- Constrast stretching with clipping (if required)
% 4- Perform a Z-projection (max or avg) (if required)

% Lighten dark regions
if lightenDarkRegions.flag
    switch lightenDarkRegions.method
        case 'log'
            cellIm = stretch(log(double(cellIm)+(10^-6)));
        case 'root'
            cellIm = stretch(double(cellIm).^lightenDarkRegions.expVal);
    end
end


% Smoothing (if required)
warning('off','MATLAB:smoothn:SLowerBound')
if isempty(smoothVal)
    [filtCell,smoothVal] = smoothn(double(cellIm));
    filtCell = dip_image(cellIm);
else
    if (smoothVal==0)
        filtCell = cellIm;
    else
        filtCell = (dip_image(smoothn(double(cellIm),smoothVal))); 
    end
end
warning('on','MATLAB:smoothn:SLowerBound')


% Contrast stretching of the cell image
filtCell = stretch(filtCell,lb,ub,0,255);
  
% Check if we have an empty image (= no cell)
try
    tmpTh= threshold(filtCell,'otsu'); 
    filtCell = im2mat(round(filtCell),'uint8');
catch % if the image can not be thresholded, it means that the image is empty
    s = size(im2mat(filtCell));
    filtCell = zeros(s,'uint8');
end



