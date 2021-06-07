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
function writeAbaqusResultsScript(write_options,coordinates,elem0,elem,u_fwd,u_inv,mech_props)


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

%Build a file name
if isfile([write_options.write_path filesep write_options.write_filename '_check.inp'])
    delete([write_options.write_path filesep write_options.write_filename '_check.inp']);
end
fid = fopen([write_options.write_path filesep write_options.write_filename '_check.inp'],'wt');

%%%%%%%%%%%%%%%%%%%%%%%% Control Abaqus inverse %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Write the heading
fprintf(fid,'*HEADING');
fprintf(fid,'\n','');

%Write the nodes
fprintf(fid,'*NODE,NSET=NTODOS');
fprintf(fid,'\n','');
nnode = 1:size(coordinates,1);
rd=[nnode', coordinates(:,1), coordinates(:,2), coordinates(:,3)];
fprintf(fid,'%7i, %12.7f, %12.7f, %12.7f\n',rd');
clear rd;

%Write the gel elements
fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=GEL');
fprintf(fid,'\n','');
rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');
clear rr;

%Write the cell elements
fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=SPROUT');
fprintf(fid,'\n','');
rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');
clear rr;

%Define:   Hydrogel = MATPROPS; cell = MATPROPS2
fprintf(fid,'*SOLID SECTION, ELSET=SPROUT, MATERIAL=MATPROPS2');
fprintf(fid,'\n','');
fprintf(fid,'*SOLID SECTION, ELSET=GEL, MATERIAL=MATPROPS');
fprintf(fid,'\n','');

%Define cell's material properties
fprintf(fid,'*MATERIAL,NAME=MATPROPS2'); 
fprintf(fid,'\n','');
fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
fprintf(fid,'\n','');
fprintf(fid,'100.0e-6,  4.50000E-01');
fprintf(fid,'\n','');

%Define hydrogel's material properties
fprintf(fid,'*MATERIAL,NAME=MATPROPS'); 
fprintf(fid,'\n','');

%Check the behavior
switch ecm_behavior    
    case 'Linear elastic' 
        %Define it as a linear elastic model
        fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
        fprintf(fid,'\n','');
        
        %Define the constants     
        % fprintf(fid,'100.0,  4.50000E-01');
        fprintf(fid,[elastic_modulus ', ' poissons_ratio]);
        fprintf(fid,'\n','');
        
        %Others
        fprintf(fid,'*STEP,INC=1000');
        fprintf(fid,'\n','');
        fprintf(fid,'*STATIC,DIRECT');
        fprintf(fid,'\n','');
        fprintf(fid,'1.0,1.0');
        fprintf(fid,'\n','');
        
    case 'Non linear elastic'
        switch constitutive_model
            case 'Neo-hookean'
                %Define it as a neo hookean model
                fprintf(fid,'*Hyperelastic, neo hooke');
                fprintf(fid,'\n','');
                
                %Define the constants
                fprintf(fid,[C10 ', ' D1]);
                fprintf(fid,'\n','');
                
                %Others
                fprintf(fid,'*STEP,NLGEOM,INC=1000');
                fprintf(fid,'\n','');
                fprintf(fid,'*STATIC');
                fprintf(fid,'\n','');
                fprintf(fid,'1.0,1.0,0.05,1.0');
                fprintf(fid,'\n','');

        end
end

indi = nnode;

abb = 1:size(coordinates,1);

%Write the forward displacement field                
fprintf(fid,'*BOUNDARY,TYPE=DISPLACEMENT');
fprintf(fid,'\n','');

rr=[abb',ones(length(indi),2),u_fwd(3*indi-2)'; abb',2*ones(length(indi),2),u_fwd(3*indi-1)'; abb', 3*ones(length(indi),2), u_fwd(3*indi)'];
fprintf(fid,'%7i, %7i, %7i, %12.7f\n',rr');
clear rr;
fprintf(fid,'*NODE FILE');
fprintf(fid,'\n','');
fprintf(fid,'U');
fprintf(fid,'\n','');
fprintf(fid,'*EL FILE,ELSET=GEL,POSITION=AVERAGED AT NODES');
fprintf(fid,'\n','');
fprintf(fid,'S,SP,EP');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');


switch ecm_behavior
    case 'Linear elastic'
        fprintf(fid,'*STEP,INC=1000');
        fprintf(fid,'\n','');
        fprintf(fid,'*STATIC,DIRECT');
        fprintf(fid,'\n','');
        fprintf(fid,'1.0,1.0');
        
    case 'Non linear elastic'
        fprintf(fid,'*STEP,NLGEOM,INC=1000');
        fprintf(fid,'\n','');
        fprintf(fid,'*STATIC');
        fprintf(fid,'\n','');
        fprintf(fid,'1.0,1.0,0.05,1.0');
end
fprintf(fid,'\n','');


%Write the inverse displacement field 
fprintf(fid,'*BOUNDARY,TYPE=DISPLACEMENT');
fprintf(fid,'\n','');

rr=[abb',ones(length(indi),2),u_inv(3*indi-2)'; abb',2*ones(length(indi),2),u_inv(3*indi-1)'; abb', 3*ones(length(indi),2), u_inv(3*indi)'];
fprintf(fid,'%7i, %7i, %7i, %12.7f\n',rr');
clear rr;

fprintf(fid,'*NODE FILE');
fprintf(fid,'\n','');
fprintf(fid,'U');
fprintf(fid,'\n','');
fprintf(fid,'*EL FILE,ELSET=GEL,POSITION=AVERAGED AT NODES');
fprintf(fid,'\n','');
fprintf(fid,'S,SP,EP');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');


fclose(fid);