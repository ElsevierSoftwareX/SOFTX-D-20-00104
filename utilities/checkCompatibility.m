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

function [abq1_text,abq2_text] = checkCompatibility()

abq1_text = [];
abq2_text = [];

[Result, WordPath] = system('WHERE /F /R "c:\Program Files (x86)" vcvarsall.bat');

if Result == 0
    tmp = splitlines(erase(convertCharsToStrings(WordPath),'"'));
    abq1_text = tmp{1};
    clear tmp;
else
    return;
end


[Result, WordPath] = system('WHERE /F /R "c:\Program Files (x86)" ifortvars.bat');
if Result == 0
    tmp = splitlines(erase(convertCharsToStrings(WordPath),'"'));
    abq2_text = tmp{1};
    clear tmp;
else
    return;
end

