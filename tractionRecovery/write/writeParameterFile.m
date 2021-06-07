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
function writeParameterFile(output_path,mesh_options,size_elem_gel,size_elem_cell,rect,mech_props)

fid = fopen([output_path filesep 'parameters.txt'], 'wt' );

fprintf(fid, ['SELECTED CHANNEL \n: ' mesh_options.selected_channel '\n']);

fprintf(fid, 'CROP PARAMETERS \n');
fprintf(fid, ['Crop: ' num2str(rect) '\n']);

fprintf(fid, 'MESH PARAMETERS \n');
fprintf(fid, ['hydrogel meshing: ' num2str(mesh_options.maxvol_gel) '\n']);
fprintf(fid, ['cell meshing: ' num2str(mesh_options.maxsize_cell) '\n']);
fprintf(fid, ['cell pre-smoothing: ' num2str(mesh_options.pre_smoothing) '\n']);
fprintf(fid, ['cell post-smoothing: ' num2str(mesh_options.post_smoothing) '\n']);
fprintf(fid, ['cell post-smoothing amount: ' num2str(mesh_options.post_smoothing_amount) '\n']);

fprintf(fid, 'RESULTANT ELEMENTS (for the last analyzed timepoint) \n');
fprintf(fid, ['number of hydrogel elements: ' num2str(size_elem_gel) '\n']);
fprintf(fid, ['number of cell elements: ' num2str(size_elem_cell,1) '\n']);
fprintf(fid, ['number of TOTAL elements: ' num2str(size_elem_cell+size_elem_gel) '\n']);

fprintf(fid, 'ECM PROPERTIES \n');
fprintf(fid, ['Behavior: ' mech_props.ecm_behavior '\n']);
fprintf(fid, ['Young`s modulus: ' num2str(mech_props.elastic_modulus) '\n']);
fprintf(fid, ['Poisson`s ratio: ' num2str(mech_props.poissons_ratio) '\n']);

fprintf(fid, 'SOLVER USED \n');
fprintf(fid, [mech_props.solver '\n']);

fclose(fid);