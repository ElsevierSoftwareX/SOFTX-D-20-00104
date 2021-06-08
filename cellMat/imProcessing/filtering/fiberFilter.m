%    This file is part of FFDcalc.
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

function [filtFibers,smoothVal] = fiberFilter(fiberIm,lightenDarkRegions,smoothVal,lb,ub)
% This function smooths out the fiber image and increase its contrast.


% Boost dark areas of the image
if lightenDarkRegions.flag
    switch lightenDarkRegions.method
        case 'log'
            fiberIm = stretch(log(double(fiberIm)+(10^-6)));
        case 'root'
            fiberIm = stretch(double(fiberIm).^lightenDarkRegions.expVal);
    end
end

% Smooth the fiber image (if required)
warning('off','MATLAB:smoothn:SLowerBound')
if isempty(smoothVal)
    [filtFibers,smoothVal] = smoothn(double(fiberIm));
    filtFibers = dip_image(filtFibers);
else
    if (smoothVal==0)
        filtFibers = fiberIm;
    else
        filtFibers = (dip_image(smoothn(double(fiberIm),smoothVal))); 
    end
end
warning('on','MATLAB:smoothn:SLowerBound')

% Contrast stretching
filtFibers = stretch(filtFibers,lb,ub,0,255);


filtFibers = im2mat(round(filtFibers),'uint8');

