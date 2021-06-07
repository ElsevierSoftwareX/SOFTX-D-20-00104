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
function writeVector2Paraview(fid,vector,vector_name,precision)
% Vector data 
fprintf(fid, ['\nVECTORS ', vector_name,' float\n']);
vecdataSpecs = [repmat(['%0.', precision, 'f '], 1, 3), '\n'];
vectorx = vector(1:3:size(vector,1));
vectory = vector(2:3:size(vector,1));
vectorz = vector(3:3:size(vector,1));
fprintf(fid, vecdataSpecs, [vectorx, vectory, vectorz]');

