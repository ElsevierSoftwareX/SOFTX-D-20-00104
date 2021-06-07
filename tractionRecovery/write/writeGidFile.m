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
function writeGidFile(coordinates,uf,ui,filename,bf,bi,tf,ti,ssf,srf,ssi,sri)

%Number the coordinates
id_coordinates = 1:size(coordinates,1);

%Displacement forward
fid = fopen(strcat(filename,'.flavia.res'),'wt');
fprintf(fid,'Gid Post Results File 1.0');
fprintf(fid,'\n','');
fprintf(fid,'RESULT "Displacements Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', uf(1:3:size(uf,1)), uf(2:3:size(uf,1)), uf(3:3:size(uf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

%Displacement inverse
fprintf(fid,'RESULT "Displacements Inverse" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', ui(1:3:size(uf,1)), ui(2:3:size(uf,1)), ui(3:3:size(uf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

%Volume forces forward
fprintf(fid,'RESULT "Volume Forces Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', bf(1:3:size(bf,1)), bf(2:3:size(bf,1)), bf(3:3:size(bf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

%Volume forces inverse
fprintf(fid,'RESULT "Volume Forces Inverse" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', bi(1:3:size(bf,1)), bi(2:3:size(bf,1)), bi(3:3:size(bf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

%Traction forces forward
fprintf(fid,'RESULT "Traction Forces Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', tf(1:3:size(tf,1)), tf(2:3:size(tf,1)), tf(3:3:size(tf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

%Traction forces inverse
fprintf(fid,'RESULT "Traction Forces Inverse" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', ti(1:3:size(tf,1)), ti(2:3:size(tf,1)), ti(3:3:size(tf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');


fprintf(fid,'RESULT "Principal Stress Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', ssf(:,1), ssf(:,2), ssf(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

fprintf(fid,'RESULT "Principal Stress Inverse" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', ssi(:,1), ssi(:,2), ssi(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

fprintf(fid,'RESULT "Principal Strain Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', srf(:,1), srf(:,2), srf(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

fprintf(fid,'RESULT "Principal Strain Inverse" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[id_coordinates', sri(:,1), sri(:,2), sri(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

fclose(fid);




