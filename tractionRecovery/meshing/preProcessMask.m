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

function mask_out = preProcessMask(mask_in,size_morphological_operator)


%Dilate
SE = strel('sphere',size_morphological_operator);
mask_in = imdilate(mask_in, SE);
%Fill holes
mask_in = imfill(mask_in,'holes');
%Erode
mask_in = imerode(mask_in, SE); clear SE;
%Keep the largest mask objects
conn = 3^ndims(mask_in)-1;
[mask_in, ~] = getLargestCc(mask_in,conn,4);
%Add padding
padSize = 3;
mask_in=padarray(mask_in,[padSize padSize padSize]);
%Smooth out the mask surface
[mask_in,~] = moisturize (mask_in);
%Remove the previous padding
mask_out = mask_in(padSize+1:end-padSize,padSize+1:end-padSize,padSize+1:end-padSize);