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

function res = getResolution(tiff_raw)

%Split the string
tmp = strsplit(tiff_raw,'_');
    
%Get the XY resolution
xy_res = tmp{3};
xy_res = str2double(xy_res(6:end));

%Return the resolution in um
res = [xy_res xy_res];

if length(tmp)>3
    
    %Get the Z resolution
    z_res = tmp{end};
    z_res = str2double(z_res(5:end));
    
    res = [res z_res];
end

