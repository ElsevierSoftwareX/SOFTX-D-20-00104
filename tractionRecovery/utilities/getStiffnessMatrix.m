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
function A = getStiffnessMatrix(result_folders,write_options,mech_props,elements,coordinates)

%Check the behavior
switch mech_props.ecm_behavior
    case 'Linear elastic'
        %Load stiffness matrix
        asparse = load([result_folders.intermediate_results_folder filesep write_options.write_filename '_STIF1.mtx']);
        %Save RAM by converting it to sparse
        A = spconvert(asparse);
        
    case 'Non linear elastic'
        %Load the stiffness matrix
        bsparse = load([result_folders.intermediate_results_folder filesep write_options.write_filename '_STIF2.mtx']);
        
        ind_i = bsparse(:,1)*12 + bsparse(:,2);
        ind_j = bsparse(:,3);        
        
        BSF = spconvert([ind_i ind_j bsparse(:,4)]);        
        clear bsparse;
        
        KS1 = (reshape(BSF',12*12,size(BSF,1)/12))';
        clear BSF;
        
        Ig = 3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
        for i=1:11
            Ig = [Ig,3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0]];
        end
        
        for i=1:4
            if i==1
                Jg = 3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[2,2,2,2,2,2,2,2,2,2,2,2];
            else
                Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[2,2,2,2,2,2,2,2,2,2,2,2]];
            end
            Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[1,1,1,1,1,1,1,1,1,1,1,1]];
            Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[0,0,0,0,0,0,0,0,0,0,0,0]];
        end
        clear elements;
        %Save RAM by converting it to sparse        
        AA = sparse ( Ig , Jg , KS1, 3*size(coordinates,1) , 3*size(coordinates,1) );
        clear Ig Jg KS1 coordinates;
        
        [n,m] = size(AA);
        A = AA' + AA;
        A(1:n+1:end) = diag(AA);
        clear AA;        
end
