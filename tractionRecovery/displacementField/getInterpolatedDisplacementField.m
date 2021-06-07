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
function [UX,UY,UZ] = getInterpolatedDisplacementField(dispField,Lx,Ly,Lz,node)

field_size = size(dispField.X);

[x,y,z] = meshgrid((1:field_size(2))-1,(1:field_size(1))-1,(1:field_size(3))-1);
x = x*Lx;
y = y*Ly;
z = z*Lz;

FUX = griddedInterpolant(y,x,z,dispField.X,'cubic'); % Build the interpolation function
UX = FUX(node(:,2),node(:,1),node(:,3));

FUY = griddedInterpolant(y,x,z,dispField.Y,'cubic'); % Build the interpolation function
UY = FUY(node(:,2),node(:,1),node(:,3));

FUZ = griddedInterpolant(y,x,z,dispField.Z,'cubic'); % Build the interpolation function
UZ = FUZ(node(:,2),node(:,1),node(:,3));


%% Old

% %Generate a grid of coordinates in um
% [rows,columns,slices] = ndgrid((1:size(dispField.X,1))*Lx,(1:size(dispField.X,2))*Ly,(1:size(dispField.X,3))*Lz);
% 
% %Create an interpolant object for X 
% % FUX = griddedInterpolant(rows,columns,slices,dispField.Y);
% FUX = griddedInterpolant(rows,columns,slices,dispField.X);
% %Interpolate the X component of the displacement field to the nodes
% UX = FUX(node(:,1),node(:,2),node(:,3));
% clear FUX
% 
% %Create an interpolant object for X 
% % FUY = griddedInterpolant(rows,columns,slices,dispField.X);
% FUY = griddedInterpolant(rows,columns,slices,dispField.Y);
% %Interpolate the X component of the displacement field to the nodes
% UY = FUY(node(:,1),node(:,2),node(:,3));
% clear FUY
% 
% %Create an interpolant object for X 
% FUZ = griddedInterpolant(rows,columns,slices,dispField.Z);
% %Interpolate the X component of the displacement field to the nodes
% UZ = FUZ(node(:,1),node(:,2),node(:,3));
% clear FUZ

