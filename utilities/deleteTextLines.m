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


function deleteTextLines(strIn,fileIn,fileOut)


if isempty(strIn{1})
    
    copyfile(fileIn,fileOut) ;
    
else
    
    fidIn = fopen(fileIn,'rt');
    fidOut = fopen(fileOut,'wt');
    while ~feof(fidIn) ;
        tline=fgets(fidIn) ;
        flag = zeros(1,length(strIn));
        % Chek if all the strings provided are missing in the current extracted line
        for ii=1:length(strIn)
            if isempty(strfind(tline, strIn{ii}))
                flag(ii) = 1;
            end
        end
        if sum(flag)==length(strIn) % if all the strings are missing, we write the current line to fileOut
            fprintf(fidOut,'%c',tline) ;
        end
        clear flag
    end
    fclose(fidIn);
    fclose(fidOut);
    
end




