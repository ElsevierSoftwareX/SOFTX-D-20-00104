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

function surfVectorField = vectorProjSurf(vectorField,surf_vert,res)

%Compute the field size
field_size = size(vectorField.X);

%Get the space coordinates
[x,y,z] = meshgrid(1:field_size(2),1:field_size(1),1:field_size(3));
x = x*res(1);
y = y*res(2);
z = z*res(3);

% Build the interpolation function
interpFun = griddedInterpolant(y,x,z,vectorField.X,'linear'); 

%Interpolate values to the surface nodes
surfVectorField.X = interpFun(surf_vert(:,2),surf_vert(:,1),surf_vert(:,3));
interpFun.Values = vectorField.Y;
surfVectorField.Y = interpFun(surf_vert(:,2),surf_vert(:,1),surf_vert(:,3));
interpFun.Values = vectorField.Z;
surfVectorField.Z = interpFun(surf_vert(:,2),surf_vert(:,1),surf_vert(:,3));




    
    