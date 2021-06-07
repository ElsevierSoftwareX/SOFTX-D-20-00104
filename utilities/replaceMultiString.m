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



function replaceMultiString(strIn,strOut,fileIn,fileOut)
% This function replaces the strings provided in the struct strIn with the
% ones in strOut. The new result is written in the file fileOut

if isempty(strIn{1})
    
    copyfile(fileIn,fileOut) ;
    
else
    
    fid = fopen(fileIn,'rt');
    a=fscanf(fid,'%c'); % Read in the whole file
    fclose(fid);
    
    for ii=1:length(strIn)
        a=strrep(a,strIn{ii},strOut{ii}); % Replace strings
    end
    
    fid = fopen(fileOut,'wt');
    fprintf(fid,'%c',a); % Write out the whole file
    fclose(fid);
    
end
