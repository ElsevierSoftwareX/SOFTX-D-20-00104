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
function writeParaviewFile(faces,coordinates,uf,ui,cell_folder,mech_vars_folder,file_name,bf,bi,tf,ti,ssf,srf,ssi,sri)

%% Cell data

 
% VTK options
precision = '2'; % of the floating point numbers
fileTitle = 'VTK from Matlab - Cell Surf';
posNum = size(coordinates,1);
faceNum = size(faces,1);

fid = fopen([cell_folder filesep 'cell_' file_name '.vtk'],'wt');
fprintf(fid, '# vtk DataFile Version 3.0\n');

% Title
fprintf(fid, [fileTitle '\n']);

% File formar (ASCII or BINARY)
fprintf(fid, 'ASCII\n');

% Define the dataset structure
fprintf(fid, 'DATASET POLYDATA\n');

%Get the faces of the mask only
faces1=faces((faces(:,4)==0),1:3);
faceNum = size(faces1,1);

% Define points (coordinates)
coordSpecs = [repmat('%0.7f ', 1, 3), '\n'];
fprintf(fid, ['POINTS ' num2str(posNum) ' float\n']);
fprintf(fid, coordSpecs, coordinates');

% Define the polygons (faces of the triangulation)
triSpecs = '%d %d %d\n';
fprintf(fid,'\nPOLYGONS %d %d\n',faceNum,4*faceNum);
fprintf(fid,['3 ' triSpecs],(faces1-1)');

% Start with point data
fprintf(fid, ['\nPOINT_DATA ' num2str(posNum)]);

% Scalar data 1 (displacement magnitude)
fprintf(fid, ['\nSCALARS ', 'dispMagFWD',' float 1\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
scaldataSpecs = [repmat(['%0.', precision, 'f '], 1, 1), '\n'];
mag = sqrt( uf(1:3:size(uf,1)).^2 + uf(2:3:size(uf,1)).^2 + uf(3:3:size(uf,1)).^2 );
fprintf(fid, scaldataSpecs, mag');

% Scalar data 1 (displacement magnitude)
fprintf(fid, ['\nSCALARS ', 'dispMagINV',' float 1\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
scaldataSpecs = [repmat(['%0.', precision, 'f '], 1, 1), '\n'];
mag = sqrt( ui(1:3:size(ui,1)).^2 + ui(2:3:size(ui,1)).^2 + ui(3:3:size(ui,1)).^2 );
fprintf(fid, scaldataSpecs, mag');

% Scalar data 1 (traction magnitude)
fprintf(fid, ['\nSCALARS ', 'tracMagFWD',' float 1\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
scaldataSpecs = [repmat(['%0.', precision, 'f '], 1, 1), '\n'];
mag = sqrt( tf(1:3:size(tf,1)).^2 + tf(2:3:size(tf,1)).^2 + tf(3:3:size(tf,1)).^2 );
fprintf(fid, scaldataSpecs, mag');

% Scalar data 1 (traction magnitude)
fprintf(fid, ['\nSCALARS ', 'tracMagINV',' float 1\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
scaldataSpecs = [repmat(['%0.', precision, 'f '], 1, 1), '\n'];
mag = sqrt( ti(1:3:size(ti,1)).^2 + ti(2:3:size(ti,1)).^2 + ti(3:3:size(ti,1)).^2 );
fprintf(fid, scaldataSpecs, mag');

fclose(fid);


clear fid;

%% Write the rest of the data 
posNum = size(coordinates,1);
fileTitle = 'VTK from Matlab - Displacement field';

fid = fopen([mech_vars_folder filesep 'mech_vars_' file_name,'.vtk'],'wt');
fprintf(fid, '# vtk DataFile Version 3.0\n');

% Title
fprintf(fid, [fileTitle '\n']);

% File formar (ASCII or BINARY)
fprintf(fid, 'ASCII\n');

% Define the dataset structure
fprintf(fid, 'DATASET POLYDATA\n');


% Define points (coordinates)
coordSpecs = [repmat('%0.7f ', 1, 3), '\n'];
fprintf(fid, ['POINTS ' num2str(posNum) ' float\n']);
fprintf(fid, coordSpecs, coordinates');



%Define a list of names of the variables that we want to write into the
%file
variable_names = {'uf','ui','tf','ti','bf','bi'};
title_names = {'dispFieldFWD','dispFieldINV','tracFieldFWD','tracFieldINV','volForFWD'...
    'volForINV'};
% Start with point data
fprintf(fid, ['\nPOINT_DATA ' num2str(posNum)]);

for ii = 1:length(variable_names)
    eval(['writeVector2Paraview(fid,',variable_names{ii},',title_names{ii},precision)']);
end

% Close the file
fclose(fid);
