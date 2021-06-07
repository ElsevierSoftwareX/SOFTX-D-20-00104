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

function experiment_name = getExperimentName(data_path)
tmp = strsplit(data_path,filesep);
tmp = strsplit(tmp{end},'_');
tmp2 = cellfun(@(x) strcmp(x,'difAnalysis'),tmp);
tmp = tmp(not(tmp2));
experiment_name = '';
for ii = 1:length(tmp)
    if ii==1
        experiment_name = [experiment_name tmp{ii}];
    else
        experiment_name = [experiment_name '_' tmp{ii}];
    end
end