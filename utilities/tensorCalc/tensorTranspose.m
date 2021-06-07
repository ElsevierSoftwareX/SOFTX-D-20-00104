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


function output = tensorTranspose(input)

dim = max(size(input));

% initialization
output = cell(size(input));
[output{:}]=deal(zeros(size(input{1})));

switch dim
    
    case 2
        if size(input,1)==size(input,2)        
            output{1,1} = input{1,1};
            output{1,2} = input{2,1};
            output{2,1} = input{1,2};
            output{2,2} = input{2,2};
        else
            output = input;
        end
        
    case 3
        if size(input,1)==size(input,2)
            output{1,1} = input{1,1};
            output{1,2} = input{2,1};
            output{1,3} = input{3,1};
            output{2,1} = input{1,2};
            output{2,2} = input{2,2};
            output{2,3} = input{3,2};
            output{3,1} = input{1,3};
            output{3,2} = input{2,3};
            output{3,3} = input{3,3};
        else
            output = input;
        end
        
end
