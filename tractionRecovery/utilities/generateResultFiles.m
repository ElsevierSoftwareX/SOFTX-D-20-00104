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

function generateResultFiles(result_folders,file_name,coordinates,u_fwd,u_inv,fc1,mechanical_variables,visualization_tool,elem_gel,elem_cell)

vc = 1:size(coordinates,1);

uf = mechanical_variables.uf;
ui = mechanical_variables.ui;
bf = mechanical_variables.bf;
bi = mechanical_variables.bi;
tf = mechanical_variables.tf;
ti = mechanical_variables.ti;
ssf = mechanical_variables.ssf;
ssi = mechanical_variables.ssi;
srf = mechanical_variables.srf;
sri = mechanical_variables.sri;

%% Write the raw data to data files

%Define the names of the variables that we will treat
variable_names = {'uf','ui','bf','bi','tf','ti','ssf','ssi','srf','sri'};
%Since vector variables have the values arranged in a different way than
%scalar variables, we need flags to know if they are vectors. There are 6
%vectors and 4 scalars.
vector_flags = [ones(1,6) zeros(1,4)]; 
vc = vc';

for jj = 1:length(variable_names)
    %Extract the variable name
    variable_name = variable_names{jj};
    %Create a result file with that name
    fid = fopen([result_folders.raw_results_folder filesep file_name '_' variable_name],'wt');    
    %Run the command to generate a variable called rd
    if vector_flags(jj)
        eval(['rd=[vc, ',variable_name,'(1:3:size(',variable_name,',1)), ',variable_name,'(2:3:size(',variable_name,',1)), ',variable_name,'(3:3:size(',variable_name,',1))];']);
    else        
        eval(['rd=[vc, ',variable_name,'(:,1), ',variable_name,'(:,2), ',variable_name,'(:,3)];']);
    end
    %     rd=[vc, uf(1:3:size(uf,1)), uf(2:3:size(uf,1)), uf(3:3:size(uf,1))];
    
    %Print the results in the file
    fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
    %Close the file
    fclose(fid);    
end

%Store the coordinates
fid = fopen([result_folders.raw_results_folder filesep file_name '_coordinates'],'wt');
rd=[vc, coordinates];
%Print the results in the file
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
%Close the file
fclose(fid);

%Store the elements
vc = 1:size(elem_gel,1);
fid = fopen([result_folders.raw_results_folder filesep file_name '_elemgel'],'wt');
rd=[vc', elem_gel];
%Print the results in the file
fprintf(fid,'%7i %12.7f %12.7f %12.7f %12.7f\n',rd');
%Close the file
fclose(fid);

%Store the elements
vc = 1:(size(elem_cell,1)+size(elem_gel,1));
fid = fopen([result_folders.raw_results_folder filesep file_name '_elements'],'wt');
rd=[elem_cell;elem_gel];
%Print the results in the file
fprintf(fid,'%7i %7i %7i %7i %7i\n',rd');
%Close the file
fclose(fid);

%% Check the chosen visualization tool

%Check which visualization tools must be selected
flags.gid = or(strcmp(visualization_tool,'gid'), strcmp(visualization_tool,'both'));
flags.paraview = or(strcmp(visualization_tool,'paraview'), strcmp(visualization_tool,'both'));

%GID
if flags.gid
    writeGidFile(coordinates,uf,ui,[result_folders.gid_results_folder filesep file_name],bf,bi,tf,ti,ssf,srf,ssi,sri);
    %Handle the case in which there is some output requested. Here the
    %sampling will be empty
    if nargout
        varargout{1} = [];
    end
end

%Paraview
if flags.paraview
        
        %Generate the folders for the paraview outputs
%         cell_folder = [result_folders.paraview_results_folder filesep 'cell'];
%         if not(exist(cell_folder,'dir')==7)
%             mkdir(cell_folder);
%         end
%         mech_vars_folder = [result_folders.paraview_results_folder filesep 'mech_vars'];
%         if not(exist(mech_vars_folder,'dir')==7)
%             mkdir(mech_vars_folder);
%         end

        if not(exist(result_folders.paraview_results_folder,'dir')==7)
            mkdir(result_folders.paraview_results_folder);
        end

        writeParaviewFile(fc1,coordinates,uf,ui,result_folders.paraview_results_folder,result_folders.paraview_results_folder,file_name,bf,bi,tf,ti,ssf,srf,ssi,sri);        
 end

