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

% assuimg that the displacement field is in the output variable:
% rows: nodes, 3 columns: u1, u2, u3
% output(:, 1) = u1;
% output(:, 2) = u2;
% output(:, 3) = u3;

function writeDisplacementFF(write_path,output,extension)

for ii = 1:3
  filename =  [write_path filesep 'u' num2str(ii) '_' extension '.sol'];
  fileID = fopen(filename, 'w');

  fprintf(fileID, 'MeshVersionFormatted 1\n\n');
  fprintf(fileID, 'Dimension 3\n\n');
  fprintf(fileID, 'SolAtVertices\n');
  fprintf(fileID, '%d\n', size(output, 1));
  fprintf(fileID, '1 1\n');

  fclose(fileID);

  dlmwrite(filename, output(:, ii), '-append', 'precision', 8);

  fileID = fopen(filename, 'a');
  fprintf(fileID, '\nEnd\n');
  fclose(fileID);
end
