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

function cell2vtk(mask,dispField,res,save_file)

%Make sure that the mask is in binary format
mask = logical(mask);

% Create a smooth surface of the cell
[~,cell_data] = moisturize(logical(mask),res);

% Project the displacements onto the surface of the cell
surf_vec_field = vectorProjSurf(dispField,cell_data.vertices,res);

%Compute the magnitude
if length(size(mask))==3
    vec_mag = sqrt(surf_vec_field.X(:).^2 + surf_vec_field.Y(:).^2 + surf_vec_field.Z(:).^2);
else
    vec_mag = sqrt(surf_vec_field.X(:).^2 + surf_vec_field.Y(:).^2);
end
%% Export to VTK

% VTK options
precision = '2'; % of the floating point numbers
fileTitle = 'VTK from Matlab - Cell Surf';
dataTitle.scalar1 =  'dispMag';

posNum = size(cell_data.vertices,1);
faceNum = size(cell_data.faces,1);

% Open the file
fid = fopen(save_file, 'w'); 

% VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 3.0\n');

% Title
fprintf(fid, [fileTitle '\n']);

% File formar (ASCII or BINARY)
fprintf(fid, 'ASCII\n');

% Define the dataset structure
fprintf(fid, 'DATASET POLYDATA\n');

% Define points (coordinates)
coordSpecs = [repmat(['%0.', precision, 'f '], 1, length(size(mask))), '\n'];
fprintf(fid, ['POINTS ' num2str(posNum) ' float\n']);
fprintf(fid, coordSpecs, cell_data.vertices');

% Define the polygons (faces of the triangulation)
triSpecs = '%d %d %d\n';
fprintf(fid,'\nPOLYGONS %d %d\n',faceNum,4*faceNum);
fprintf(fid,['3 ' triSpecs],(cell_data.faces-1)');

% Start with point data
fprintf(fid, ['\nPOINT_DATA ' num2str(posNum)]);

% Scalar data 1 (displacement magnitude)
fprintf(fid, ['\nSCALARS ', dataTitle.scalar1,' float 1\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
scaldataSpecs = [repmat(['%0.', precision, 'f '], 1, 1), '\n'];
fprintf(fid, scaldataSpecs, vec_mag');

% Close the file
fclose(fid);

