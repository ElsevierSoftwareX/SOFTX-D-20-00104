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


function im = shiftCorr(im,shft,cropStart,cropEnd,type)


% Change the boundary condtions to padd with zeros the locations where the
% image won't exist due to the shift correction
oldBoundaryOption = dip_getboundary(1);
dip_setboundary('add_zeros');


% Apply the absolute shift correction
shft = -shft; % the sign has to be inverted to be used with the "resample" function
switch type
    case 'bin'
        im = resample(dip_image(im)>0,1,shft,'nn');
    case {'field','im'}
        im = resample(dip_image(im),1,shft,'bspline');
end

% Crop
imDim = length(size(im));
switch imDim
    case 2
        im = im(cropStart(1):end-cropEnd(1),cropStart(2):end-cropEnd(2));
    case 3
        im = im(cropStart(1):end-cropEnd(1),cropStart(2):end-cropEnd(2),cropStart(3):end-cropEnd(3));
end

switch type
    case 'bin'
        im = uint8(im2mat(1*im));
    case 'im'
        im = uint8(im2mat(round(stretch(im))));
    case 'field'
        im = double(im);
end


% Change the boundary conditions back to its previous value
dip_setboundary(oldBoundaryOption);
