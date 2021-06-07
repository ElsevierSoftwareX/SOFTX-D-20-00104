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


function output = tensorDeterminant(input)

% It is assumed that input is a square tensor of size 2x2 or 3x3


dim = max(size(input));


switch dim
    
    case 2
        
        output = (input{1,1}.*input{2,2}) - (input{1,2}.*input{2,1});

        
    case 3
        
        d1 = input{1,1}.*((input{3,3}.*input{2,2}) - (input{3,2}.*input{2,3}));
        d2 = input{2,1}.*((input{3,3}.*input{1,2}) - (input{3,2}.*input{1,3}));
        d3 = input{3,1}.*((input{2,3}.*input{1,2}) - (input{2,2}.*input{1,3}));
        output = d1 - d2 + d3;
 
end
