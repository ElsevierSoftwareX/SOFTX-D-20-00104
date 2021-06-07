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
function  diffCell2vtk(diffcellmask,resolution,diffParam,saveFile)


%% Get the data

data = diffSurfCalc(diffcellmask,resolution);
colorMtx = zeros(2,4);
colorMtx(1,:) = [diffParam.color.prot/255, 1.0];
colorMtx(2,:) = [diffParam.color.retrac/255, 1.0];



%% Export to VTK


% VTK options
precision = '2'; % of the floating point numbers
fileTitle = 'VTK from Matlab - Diff Cell';
dataTitle.scalar =  'diffCell';
posNum = size(data.vertices,1);
faceNum = size(data.faces,1);



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
fprintf(fid, coordSpecs, data.vertices');

% Define the polygons (faces of the triangulation)
triSpecs = '%d %d %d\n';
fprintf(fid,'\nPOLYGONS %d %d\n',faceNum,4*faceNum);
fprintf(fid,['3 ' triSpecs],(data.faces-1)');


% Start with point data
fprintf(fid, ['\nPOINT_DATA ' num2str(posNum)]);

% Scalar data (color data)
fprintf(fid, ['\nSCALARS ', dataTitle.scalar,' int 1\n']);
fprintf(fid, '\nLOOKUP_TABLE colorTable\n');
scaldataSpecs = [repmat('%d ', 1, 1), '\n'];
fprintf(fid, scaldataSpecs, data.colorIndx');
fprintf(fid, '\nLOOKUP_TABLE colorTable 2\n'); % define the colors for each color indx
colorSpecs = [repmat(['%0.', precision, 'f '], 1, 4), '\n'];
fprintf(fid, colorSpecs, colorMtx');

% Close the file
fclose(fid);




