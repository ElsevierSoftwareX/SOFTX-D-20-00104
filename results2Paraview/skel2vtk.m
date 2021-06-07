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

function skel2vtk(skel,mask,resolution,saveFile)


%% Get the data

% Triangulation of the cell surface and corresponding normals
th =0.1;
filter_flag = true;

%Dilate a bit the skeleton
skel = logical(im2mat(bdilation(skel,2)));
skel = skel & mask;
% skel = im2mat(gaussf(skel,1));
[skelData,~] = getSurfNormals(skel,th,resolution,filter_flag);

%% Export to VTK

% VTK options
precision = '2'; % of the floating point numbers
fileTitle = 'VTK from Matlab - Skel Surf';

posNum = size(skelData.vertices,1);
faceNum = size(skelData.faces,1);

% Open the file
fid = fopen(saveFile, 'w'); 

% VTK DataFile Version
fprintf(fid, '# vtk DataFile Version 3.0\n');

% Title
fprintf(fid, [fileTitle '\n']);

% File formar (ASCII or BINARY)
fprintf(fid, 'ASCII\n');

% Define the dataset structure
fprintf(fid, 'DATASET POLYDATA\n');

% Define points (coordinates)
coordSpecs = [repmat(['%0.', precision, 'f '], 1, 3), '\n'];
fprintf(fid, ['POINTS ' num2str(posNum) ' float\n']);
fprintf(fid, coordSpecs, skelData.vertices');

% Define the polygons (faces of the triangulation)
triSpecs = '%d %d %d\n';
fprintf(fid,'\nPOLYGONS %d %d\n',faceNum,4*faceNum);
fprintf(fid,['3 ' triSpecs],(skelData.faces-1)');

% Start with point data
fprintf(fid, ['\nPOINT_DATA ' num2str(posNum)]);

% Close the file
fclose(fid);

