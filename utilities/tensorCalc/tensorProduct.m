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


function output = tensorProduct(input1,input2)

% Important. The order of the product is : input1 x input2



% input1 and input2 can be 2 different types of tensors:
% I. A square tensor of size 2x2 or 3x3 
% II. An array-like tensor of size 1x2 (or 2x1) or 1x3 (or 3x1)

dim = max(size(input2));


% Compute  the sizes of the tensor input1 and input2
if size(input1,1)==size(input1,2)
    if size(input2,1)==size(input2,2)
        tensorDim =0;
    else
         tensorDim =1;
    end
else
   if size(input2,1)==size(input2,2)
       tensorDim =2;
   else
        tensorDim =3;
   end
end

switch dim
    
    case 2
        
        switch tensorDim     
            case 0 % input1 and input2 are square tensors of size 2x2
                output{1,1} = (input1{1,1}.*input2{1,1}) + (input1{1,2}.*input2{2,1});
                output{1,2} = (input1{1,1}.*input2{1,2}) + (input1{1,2}.*input2{2,2});
                output{2,1} = (input1{2,1}.*input2{1,1}) + (input1{2,2}.*input2{2,1});
                output{2,2} = (input1{2,1}.*input2{1,2}) + (input1{2,2}.*input2{2,2});          
            case 1 % input1 is a square tensor of size 2x2 and input2 is an array-like tensor of size 2x1 (input2 is a column tensor)
                output{1,1} = (input1{1,1}.*input2{1}) + (input1{1,2}.*input2{2});
                output{2,1} = (input1{2,1}.*input2{1}) + (input1{2,2}.*input2{2});
            case 2 % input1 is an array-like tensor of size 1x2 and input2 is a sqare tensor of siz 2x2 (input1 is a row tensor)
                output{1,1} = (input1{1}.*input2{1,1}) + (input1{2}.*input2{2,1});
                output{1,2} = (input1{1}.*input2{1,2}) + (input1{2}.*input2{2,2});
            case 3 % input1 and input 2 are array-like tensors of size 1x2 (2x1) (input1 is column tensor and input2 is a row tensor)
                output{1,1} = (input1{1}.*input2{1}) ;
                output{1,2} = (input1{1}.*input2{2}) ;
                output{2,1} = (input1{2}.*input2{1}) ;
                output{2,2} = (input1{2}.*input2{2}) ;
        end

    case 3
        
        switch tensorDim
            case 0 % input1 and input2 are square tensors of size 3x3
                output{1,1} = (input1{1,1}.*input2{1,1}) + (input1{1,2}.*input2{2,1}) + (input1{1,3}.*input2{3,1});
                output{1,2} = (input1{1,1}.*input2{1,2}) + (input1{1,2}.*input2{2,2}) + (input1{1,3}.*input2{3,2});
                output{1,3} = (input1{1,1}.*input2{1,3}) + (input1{1,2}.*input2{2,3}) + (input1{1,3}.*input2{3,3});
                output{2,1} = (input1{2,1}.*input2{1,1}) + (input1{2,2}.*input2{2,1}) + (input1{2,3}.*input2{3,1});
                output{2,2} = (input1{2,1}.*input2{1,2}) + (input1{2,2}.*input2{2,2}) + (input1{2,3}.*input2{3,2});
                output{2,3} = (input1{2,1}.*input2{1,3}) + (input1{2,2}.*input2{2,3}) + (input1{2,3}.*input2{3,3});
                output{3,1} = (input1{3,1}.*input2{1,1}) + (input1{3,2}.*input2{2,1}) + (input1{3,3}.*input2{3,1});
                output{3,2} = (input1{3,1}.*input2{1,2}) + (input1{3,2}.*input2{2,2}) + (input1{3,3}.*input2{3,2});
                output{3,3} = (input1{3,1}.*input2{1,3}) + (input1{3,2}.*input2{2,3}) + (input1{3,3}.*input2{3,3});
            case 1 %  input1 is a square tensor of size 3x3 and input2 is an array-like tensor of size 3x1 (input2 is a column tensor)
                output{1,1} = (input1{1,1}.*input2{1}) + (input1{1,2}.*input2{2}) + (input1{1,3}.*input2{3});
                output{2,1} = (input1{2,1}.*input2{1}) + (input1{2,2}.*input2{2}) + (input1{2,3}.*input2{3});
                output{3,1} = (input1{3,1}.*input2{1}) + (input1{3,2}.*input2{2}) + (input1{3,3}.*input2{3});
            case 2 % input1 is an array-like tensor of size 1x3 and input2 is a sqare tensor of siz 3x3 (input1 is a row tensor)
                output{1,1} = (input1{1}.*input2{1,1}) + (input1{2}.*input2{2,1}) + (input1{3}.*input2{3,1});
                output{1,2} = (input1{1}.*input2{1,2}) + (input1{2}.*input2{2,2}) + (input1{3}.*input2{3,2});
                output{1,3} = (input1{1}.*input2{1,3}) + (input1{2}.*input2{2,3}) + (input1{3}.*input2{3,3});
            case 3 % input1 and input 2 are array-like tensors of size 1x3 (3x1) (input1 is column tensor and input2 is a row tensor)
                output{1,1} = (input1{1}.*input2{1}); 
                output{1,2} = (input1{1}.*input2{2}); 
                output{1,3} = (input1{1}.*input2{3}); 
                output{2,1} = (input1{2}.*input2{1}); 
                output{2,2} = (input1{2}.*input2{2});
                output{2,3} = (input1{2}.*input2{3});
                output{3,1} = (input1{3}.*input2{1}); 
                output{3,2} = (input1{3}.*input2{2});
                output{3,3} = (input1{3}.*input2{3});
        end
        
end