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

function disp2vtk(dispField,resolution,sampling_indx,channel_name,save_file)

%Compute the field size
field_size = size(dispField.X);

%Get the coordinates
coords = getCoordsForVTK(field_size,resolution,sampling_indx);

% Keep the coordinates and the data at requested sampling locations
if length(field_size) == 3    
    data.vector = [dispField.X(sampling_indx(:)),dispField.Y(sampling_indx(:)),dispField.Z(sampling_indx(:))];
else
    data.vector = [dispField.X(sampling_indx(:)),dispField.Y(sampling_indx(:))];
end
clear dispField;
%% Export to VTK

% VTK options
precision = '2'; % of the floating point numbers
fileTitle = 'VTK from Matlab - Displacement Field';
dataTitle.vector =  ['dispField_' channel_name];
posNum = size(coords,1);

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
coordSpecs = [repmat(['%0.', precision, 'f '], 1, length(field_size)), '\n'];

fprintf(fid, ['POINTS ' num2str(posNum) ' float\n']);

fprintf(fid, coordSpecs, coords');

% Start with point data
fprintf(fid, ['\nPOINT_DATA ' num2str(posNum)]);

% Vector data (displacements)
fprintf(fid, ['\nVECTORS ', dataTitle.vector,' float\n']);
vecdataSpecs = [repmat(['%0.', precision, 'f '], 1, 3), '\n'];
fprintf(fid, vecdataSpecs, data.vector');

% Close the file
fclose(fid);