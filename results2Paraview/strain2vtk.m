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

% function strain2vtk(pC,volStrain,resolution,samplingIndx,saveFile)
 function strain2vtk(vol_strain,res,sampling_indx,channel_name,save_file)

%% Build the coordinates for the whole field size
field_size = size(vol_strain);

%Get the coordinates
coords = getCoordsForVTK(field_size,res,sampling_indx);

%% Export to VTK

% VTK options
precision = '2'; % of the floating point numbers
fileTitle = 'VTK from Matlab - Strain Fields';

data_title =  ['volStrain_' channel_name];

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

% Scalar data (volStrain)
fprintf(fid, ['\nSCALARS ', data_title,' float\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
scaldataSpecs = [repmat('%d ', 1, 1), '\n'];
fprintf(fid, scaldataSpecs, vol_strain(sampling_indx(:))');


% Close the file
fclose(fid);