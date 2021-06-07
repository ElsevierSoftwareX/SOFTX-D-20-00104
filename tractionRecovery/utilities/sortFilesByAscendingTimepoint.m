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
function files = sortFilesByAscendingTimepoint(files,input_type)

%Compute the number of files
n_files = length(files);
if strcmp(input_type, 'cellMat')
    n_files = n_files-1; %We remove the tR
end

%Define a container for the names of the files
file_names = cell(n_files,1);

%Define a container for the timepoint number
time_points = zeros(n_files,1);

for ii = 1:n_files
    
    %Store the file name
    file_names{ii} = files(ii).name;    

    %Split the name of the file by _
    tmp1 = strsplit(files(ii).name,'_');
        
    %Look for the one that has the word t
    tmp2 = cellfun(@(x) strfind(x(1),'t'),tmp1,'UniformOutput',false);
    tmp2 = not(cellfun(@isempty,tmp2));
    tmp1 = tmp1{tmp2};
    clear tmp2;
    
    %Get the timepoint number
    time_points(ii) = str2double(regexp(tmp1,'\d*','Match'));
    
end
clear tmp1;

%Sort the files according to the timepoints
[~,indx]=sort(time_points);
files = file_names(indx);

