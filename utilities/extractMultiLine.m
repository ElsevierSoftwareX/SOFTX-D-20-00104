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



function lineOut = extractMultiLine(strIn,fileIn)
% This function extract the lines from a text file that contain the
% provided keywords in the struct strIn

lineOut = cell(1,length(strIn));
fidIn = fopen(fileIn,'rt');
while ~feof(fidIn) ;
    tline=fgets(fidIn) ; 
    for ii=1:length(strIn)
        flag = not(isempty(strfind(tline, strIn{ii})));
        if flag
            lineOut{ii} = tline;
        end
    end

end
fclose(fidIn);

