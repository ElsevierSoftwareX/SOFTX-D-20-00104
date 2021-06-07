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

function coords = getCoordsForVTK(field_size,resolution,sampling_indx)

% Build the coordinates for the whole field size
if length(field_size) == 3
    [x,y,z] = meshgrid(1:field_size(2),1:field_size(1),1:field_size(3));
    x = x*resolution(1);
    y = y*resolution(2);
    z = z*resolution(3);
    
    coords = [x(sampling_indx(:)),y(sampling_indx(:)),z(sampling_indx(:))];
    clear x y z;
    
    
else
    [x,y] = meshgrid(1:field_size(2),1:field_size(1));
    x = x*resolution(1);
    y = y*resolution(2);
    
    coords = [x(sampling_indx(:)),y(sampling_indx(:))];
    clear x y;
    
end