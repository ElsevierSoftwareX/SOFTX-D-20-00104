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
function writeAbaqusStMtxScript(coordinates,elem0,elem,write_options,mech_props,u_fwd)

%Get the behavior parameters
ecm_behavior = mech_props.ecm_behavior;
constitutive_model = mech_props.constitutive_model;

switch ecm_behavior
    case 'Linear elastic'
        %Get the elastic modulus and the Poisson's ratio
        elastic_modulus = num2str(mech_props.elastic_modulus);
        poissons_ratio = num2str(mech_props.poissons_ratio);
    case 'Non linear elastic'
        switch constitutive_model
            case 'Neo-hookean'
                %Calculate the temperature dependent material parameters (for Neo-Hookean)
                [C10,D1] = getNeoHookeanParameters(mech_props.elastic_modulus,mech_props.poissons_ratio);
                C10 = num2str(C10);
                D1 = num2str(D1);                
        end
end

%Define a name for the abaqus file
% abaqus_file_name = [write_options.write_filename '_void'];

%Create a .inp file
if isfile([write_options.write_path filesep write_options.write_filename '.inp'])
    delete([write_options.write_path filesep write_options.write_filename '.inp']);
end
fid = fopen([write_options.write_path filesep write_options.write_filename '.inp'],'wt');

fprintf(fid,'*HEADING');
fprintf(fid,'\n','');
fprintf(fid,'*NODE,NSET=NTODOS');
fprintf(fid,'\n','');

nnode = 1:size(coordinates,1);

rd=[nnode', coordinates(:,1), coordinates(:,2), coordinates(:,3)];
fprintf(fid,'%7i, %12.7f, %12.7f, %12.7f\n',rd');

switch mech_props.abaqus.element_type
    case 'C3D4'
        
        %Hydrogel
        fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=GEL');
        fprintf(fid,'\n','');
        rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
        fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');
        
        %Sprout
        fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=SPROUT');
        fprintf(fid,'\n','');
        clear rr;
        rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
        fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');
        
    case 'C3D10'
        %Hydrogel
        fprintf(fid,'*ELEMENT,TYPE=C3D10,ELSET=GEL');
        fprintf(fid,'\n','');
        rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4),...
            elem0(:,5), elem0(:,6), elem0(:,7), elem0(:,8),elem0(:,9), elem0(:,10)];
        fprintf(fid,'%7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i\n',rr');
        
        %Sprout
        fprintf(fid,'*ELEMENT,TYPE=C3D10,ELSET=SPROUT');
        fprintf(fid,'\n','');
        clear rr;
        rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4),...
            elem(:,5), elem(:,6), elem(:,7), elem(:,8),elem(:,9), elem(:,10)];
        fprintf(fid,'%7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i, %7i\n',rr');
end

fprintf(fid,'*SOLID SECTION, ELSET=SPROUT, MATERIAL=MATPROPS2');
fprintf(fid,'\n','');
fprintf(fid,'*SOLID SECTION, ELSET=GEL, MATERIAL=MATPROPS');
fprintf(fid,'\n','');
fprintf(fid,'*MATERIAL,NAME=MATPROPS2');
fprintf(fid,'\n','');
fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
fprintf(fid,'\n','');
fprintf(fid,'100.0e-6,  4.50000E-01');
% fprintf(fid,'100.0e-6,  4.00000E-01');
fprintf(fid,'\n','');
fprintf(fid,'*MATERIAL,NAME=MATPROPS');
fprintf(fid,'\n','');

%Check the behavior
switch ecm_behavior
    case 'Linear elastic'
        fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
        fprintf(fid,'\n','');
        % fprintf(fid,'100.0,  4.50000E-01');
        fprintf(fid,[elastic_modulus ', ' poissons_ratio]);
        fprintf(fid,'\n','');
        fprintf(fid,'*STEP');
        fprintf(fid,'\n','');
        fprintf(fid,'*MATRIX GENERATE, STIFFNESS');
        fprintf(fid,'\n','');
        fprintf(fid,'*MATRIX OUTPUT, STIFFNESS, FORMAT=COORDINATE');
        fprintf(fid,'\n','');
        fprintf(fid,'*END STEP');
        fprintf(fid,'\n','');
        
    case 'Non linear elastic'
        switch constitutive_model
            case 'Neo-hookean'
                fprintf(fid,'*Hyperelastic, neo hooke');
                fprintf(fid,'\n','');              
                fprintf(fid,[C10 ', ' D1]);
                fprintf(fid,'\n','');
                fprintf(fid,'*STEP,NLGEOM,INC=1000');
                fprintf(fid,'\n','');
                fprintf(fid,'*STATIC');
                fprintf(fid,'\n','');
                fprintf(fid,'1.0,1.0,0.025,1.0');
                fprintf(fid,'\n','');
                fprintf(fid,'*BOUNDARY,TYPE=DISPLACEMENT');
                fprintf(fid,'\n','');
                
                clear rr;
                
                indi = nnode;
                
                abb = 1:size(coordinates,1);
                
                rr=[abb',ones(length(indi),2),u_fwd(3*indi-2)'; abb',2*ones(length(indi),2),u_fwd(3*indi-1)'; abb', 3*ones(length(indi),2), u_fwd(3*indi)'];
                fprintf(fid,'%7i, %7i, %7i, %12.7f\n',rr');
                fprintf(fid,'*NODE FILE');
                fprintf(fid,'\n','');
                fprintf(fid,'RF');
                fprintf(fid,'\n','');
                fprintf(fid,'*END STEP');
                fprintf(fid,'\n','');
                fprintf(fid,'*STEP,NLGEOM');
                fprintf(fid,'\n','');
                fprintf(fid,'*MATRIX GENERATE, STIFFNESS, ELEMENT BY ELEMENT');
                fprintf(fid,'\n','');
                fprintf(fid,'*MATRIX OUTPUT, STIFFNESS, FORMAT=COORDINATE');
                fprintf(fid,'\n','');
                fprintf(fid,'*END STEP');
        end
end

fclose(fid);
