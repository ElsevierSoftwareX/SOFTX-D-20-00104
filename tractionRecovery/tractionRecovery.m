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

function tractionRecovery(store_path,folder_names,file_names,experiment_name,mech_props,mesh_options,abq1_text,abq2_text,code_names)

%% User parameters
f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Traction calculation';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');

%Define the solver that you want to use to compute the tractions.
% Option 'gmres' uses the iterative solver from Matlab. It is a slow
% algorithm but it has lower RAM demands. Option 'suitesparse' uses the
% direct solver QR from the library SuiteSparse. Much faster but it
% requires a large RAM.
solver_name = 'suitesparse'; % 'suitesparse' or 'gmres'

% Define the visualization tool for which you want the outputs
visualization_tool = 'both'; %'gid', 'paraview' or 'both'

%Get the channel that we need to use
channel_name = mesh_options.selected_channel;

%Create the traction results folder and store the parameter file
save_path = [store_path filesep experiment_name filesep folder_names.tractions filesep channel_name];
if not(isfolder(save_path))
    mkdir(save_path);
end
save([save_path filesep file_names.trac_param],'mech_props','mesh_options');

%Compute the number of time points
files = dir([store_path filesep experiment_name filesep folder_names.displacements filesep channel_name]);
files = files(3:end);
nTP = length(files);

% Define the names of the folders that contain the utilities for abaqus and
% Paraview
tools_abaqus_folder_name = ['tractionRecovery' filesep 'tools_abaqus'];
tools_freefem_folder_name = ['tractionRecovery' filesep 'tools_freefem'];
tools_paraview_folder_name = ['tractionRecovery' filesep 'tools_paraview'];

%Define the name of the utilities files
utilities_names = {'abq1.bat','abq2.bat','wrt_us.for','wrt_rf.for'};


%Define the name of the folder where different results will be stored
result_folders.intermediate_results_folder = [save_path filesep folder_names.intermediate_results_folder];
if not(exist(result_folders.intermediate_results_folder,'dir')==7)
    mkdir(result_folders.intermediate_results_folder);
end


%% Traction recovery
%Load the cropping rectangle for the selected channel 
load([store_path filesep experiment_name filesep folder_names.results_for_paraview filesep channel_name filesep file_names.paraview_param],'crop_rect');

%Get the spatial resolution
res = getResolution(folder_names.tiff_raw);

%Loop through the time points 
% nTP=1;
for tt = 1:nTP
    
    %Get the file name
    file_name = files(tt).name;
    
    %Update progress bar
    d.Value = tt/nTP;
    
    %Load the resolution and mask data
    d.Message = sprintf(['Channel ' channel_name ': loading the cell mask of ' file_name(1:end-4) '...']);
    
    %Load the mask
    load([store_path filesep experiment_name filesep folder_names.post_shift_corrected filesep channel_name filesep folder_names.cellSeg filesep file_name],'cell_seg');
    
    %Crop the cell mask (the previous use of crop_rect was used in the main
    %just for the tractionRecovery GUI!! Those cropped masks were no longer
    %used)
    cell_seg = cell_seg(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
    
    
    if mesh_options.pre_smoothing > 0       
        
        %Smooth a bit because normally is very spiky
        N = mesh_options.pre_smoothing;
        kernel = ones(N, N, N) / N^3;
        blurry_image = convn(double(cell_seg), kernel, 'same');
        
        %Make it binary again
        cell_seg = uint8(blurry_image > 0.5);
    end
    
    %     cell_seg = padarray(cell_seg,[5,5,5]);
    %     cell_seg = imgaussfilt3(cell_seg,1)>0.1;% extra smooth on
    % %     cell_seg = im2mat(gaussf(cell_seg,1))>0.1;
    %     cell_seg = cell_seg(5+1:end-5,5+1:end-5,5+1:end-5);      
    
    %Interpolate the displacement field
    d.Message = sprintf(['Interpolating ' channel_name ' displacement field to the nodes of ' file_name(1:end-4) '...']);
    
    %Load the displacement field
    load([files(tt).folder filesep file_name],'dispField');
    
    %Crop the dispacement field
    dispField.X = dispField.X(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
    dispField.Y = dispField.Y(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
    dispField.Z = dispField.Z(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
         
    %Get the nodes of the mask
    d.Message = sprintf(['Computing mask nodes of ' file_name(1:end-4) '...']);
    write_options.write_path = result_folders.intermediate_results_folder;
    write_options.write_filename = file_name(1:end-4);
        

    switch mech_props.ecm_behavior
        case 'Linear elastic'
            [coordinates,elem_gel,elem_cell,faces,mesh_options] = getMeshedVolume(cell_seg,res(1),res(2),res(3),mesh_options,write_options);
        case 'Non linear elastic'
            %Get the FE mesh
            [coordinates,elem_gel,elem_cell,faces,mesh_options] = getMeshedVolume2(cell_seg,res(1),res(2),res(3),mesh_options,mech_props,write_options,dispField);
            
            %Get the interior voxels of the cell
%             cell_seg_interior = logical(im2mat(berosion(logical(cell_seg))));
            
            %Remove the displacement field values from the voxels inside
            %the cell
%             dispField.X(cell_seg_interior) = nan;
%             dispField.Y(cell_seg_interior) = nan;
%             dispField.Z(cell_seg_interior) = nan;
            
            %Interpolate
%             dispField.X = inpaintn(dispField.X,20);
%             dispField.Y = inpaintn(dispField.Y,20);
%             dispField.Z = inpaintn(dispField.Z,20);            
            
    end
    clear cell_seg;
    
    [UX,UY,UZ] = getInterpolatedDisplacementField(dispField,res(1),res(2),res(3),coordinates);
    clear dispField;
    
    %Rearrange the displacement field into a matrix
    u_fwd = zeros(1,3*size(coordinates,1));
    u_fwd(1:3:3*size(coordinates,1)) = UX;
    u_fwd(2:3:3*size(coordinates,1)) = UY;
    u_fwd(3:3:3*size(coordinates,1)) = UZ;
    clear UX UY UZ;
    
%     if strcmp(mech_props.ecm_behavior,'Non linear elastic')        
%         %"Relax" the mesh
%         coordinates_stressed = coordinates;
%         coordinates = coordinates - [u_fwd(1:3:3*size(coordinates,1))',...
%             u_fwd(2:3:3*size(coordinates,1))',u_fwd(3:3:3*size(coordinates,1))'];
%         
%     end
%     
%      figure
%      plotmesh(coordinates,faces);
%      figure
%      plotmesh(coordinates_stressed,faces);
    
    %Define results folder names
    result_folders.abaqus_results_folder = [save_path filesep folder_names.abaqus_results_folder];
    result_folders.raw_results_folder = [save_path filesep folder_names.raw_results_folder];
    result_folders.gid_results_folder = [save_path filesep folder_names.gid_results_folder];
    result_folders.paraview_results_folder = [save_path filesep folder_names.paraview_results_folder];
    result_folders.freefem_results_folder = [save_path filesep folder_names.freefem_results_folder];
    
    %Check if the result folders exist, otherwise create them
    tmp = fieldnames(result_folders);
    for jj = 1:length(tmp)
        if not(exist(result_folders.(tmp{jj}),'dir')==7)
            mkdir(result_folders.(tmp{jj}));
        end
    end
    clear jj tmp;
        
    %Write the data to a file that can be read by Abaqus
    d.Message = sprintf(['Creating ' mech_props.solver ' data file of ' file_name(1:end-4) '...']);
    write_options.write_path = result_folders.intermediate_results_folder;
    
    switch mech_props.ecm_behavior
        case 'Linear elastic'
            switch mech_props.solver
                case 'Abaqus'
                    writeAbaqusStMtxScript(coordinates,elem_gel,elem_cell,write_options,mech_props,u_fwd);
                case 'FreeFEM'
                    %Remove the 0 labels for 2 (just for the savemesh function)
                    tmp = elem_gel;
                    tmp(elem_gel(:,5)==0,5) = 2;
                    %Write the .msh file
%                     savemsh(coordinates,[tmp;elem_cell],[write_options.write_path filesep file_name(1:end-4) '_mesh.msh']);
                    savemsh(coordinates,[elem_cell;tmp],[write_options.write_path filesep file_name(1:end-4) '_mesh.msh']);

                    clear tmp;
            end
        case 'Non linear elastic'
            new_coords = coordinates-[u_fwd(1:3:end)', u_fwd(2:3:end)', u_fwd(3:3:end)'];
            writeAbaqusStMtxScript(new_coords,elem_gel,elem_cell,write_options,mech_props,u_fwd);
    end
    
    %Obtain the stiffness matrix
    switch mech_props.solver
        case 'Abaqus'
            %Save the text lines
            write_options.abq1_text = abq1_text;
            write_options.abq2_text = abq2_text;            
            %Move the tools that we need to the results path so that they can be run
            %along with the _check file
            for jj = 1:length(utilities_names)
                copyfile([tools_abaqus_folder_name filesep utilities_names{jj}],write_options.write_path);
            end            
            %Run Abaqus
            d.Message = sprintf(['Running Abaqus to obtain stiffness matrix of ' file_name(1:end-4) '...']);
            runAbaqusFile(write_options,mech_props,'stiffness_matrix');            
            %Load stiffness matrix
            d.Message = sprintf(['Loading stiffness matrix of ' file_name(1:end-4) '...']);
            A = getStiffnessMatrix(result_folders,write_options,mech_props,[elem_gel;elem_cell],coordinates);
            
        case 'FreeFEM'
            %% Do the same with FREE FEM     
            command_line = ['FreeFem++ ' [tools_freefem_folder_name filesep 'stiffness.edp']...
                ' -mesh ' [write_options.write_path filesep file_name(1:end-4) '_mesh.msh'] ...
                ' -out ' [write_options.write_path filesep file_name(1:end-4) '_matrix.txt'] ...
                ' -ngel 2'...
                ' -ncell 1'...
                ' -egel ' num2str(mech_props.elastic_modulus) ...
                ' -nugel ' num2str(mech_props.poissons_ratio)];            
            status = system(command_line);            
            A = getStiffnessMatrixFF([write_options.write_path filesep file_name(1:end-4) '_matrix.txt']);
            %Change the sign
            A = -A;
    end
    %%
    
    %Compute the forces
    switch mech_props.ecm_behavior
        case 'Linear elastic'
            %Solve the forward problem
            d.Message = sprintf(['Solving the forward problem of ' file_name(1:end-4) '...']);
            t_fwd = A*u_fwd';
            
            %Solve the inverse problem
            d.Message = sprintf(['Solving the inverse problem of ' file_name(1:end-4) '...']);
            [u_inv,elapsed_time_inv] = calculateInverseDisplacementField(elem_cell,faces,u_fwd,A,solver_name);
            t_inv = A*u_inv';
            
        case 'Non linear elastic'
            %Solve the forward problem
            d.Message = sprintf(['Solving the forward problem of ' file_name(1:end-4) '...']);
            t0_aux = load([result_folders.intermediate_results_folder filesep 'reac_for.dat']);
            t_fwd = reshape(t0_aux(:,2:end)',3*size(coordinates,1),1);
            clear t0_aux;
            
            %Solve the inverse problem
            d.Message = sprintf(['Solving the inverse problem of ' file_name(1:end-4) '...']);
            [t_inv,u_inv,elapsed_time_inv] = calculateNonLinInverseDisplacementField(result_folders,write_options,mech_props,elem_gel,elem_cell,coordinates,faces,solver_name,u_fwd);
    end
    
    disp(['Solved the inverse problem in ' num2str(floor(elapsed_time_inv/60)) ' minutes and ' num2str(round(elapsed_time_inv - floor(elapsed_time_inv/60) * 60)) ' seconds.']);
    
    %Write the Abaqus script to generate all the mechanical problem results
    d.Message = sprintf(['Computing mechanical variables of ' file_name(1:end-4) '...']);
            
    switch mech_props.solver
        case 'Abaqus'
            %Delete the utilities
            for jj = 1:length(utilities_names)
                delete([write_options.write_path filesep utilities_names{jj}]);
            end
            
            write_options.write_path = result_folders.abaqus_results_folder;
            writeAbaqusResultsScript(write_options,coordinates,elem_gel,elem_cell,u_fwd,u_inv,mech_props);
            
            %Move the tools that we need to the results path so that they can be run
            %along with the _check file
            for jj = 1:length(utilities_names)
                copyfile([tools_abaqus_folder_name filesep utilities_names{jj}],write_options.write_path);
            end
            
            %Run the Abaqus file to generate all the mechanical problem results
            runAbaqusFile(write_options,mech_props,'results');
            
            %Delete the utilities
            for jj = 1:length(utilities_names)
                delete([write_options.write_path filesep utilities_names{jj}]);
            end
            
            %Gerate all the necessary mechanical variables
            mechanical_variables = generateMechanicalVariables(write_options.write_path,elem_gel,elem_cell,u_fwd,u_inv,t_fwd,t_inv,coordinates,faces,mech_props.solver);
            
        case 'FreeFEM'
            
            %% DO THE SAME WITH FREE FEM
            write_options.write_path = result_folders.intermediate_results_folder;
            
            %Write the displacements from the forward method
            writeDisplacementFF(write_options.write_path,[u_fwd(1:3:end)', u_fwd(2:3:end)', u_fwd(3:3:end)'],'fwd');  
            
            %Calculate the stresses
            runFFresultFile(tools_freefem_folder_name,write_options,file_name,result_folders,code_names,mech_props,'fwd');
            
            %Write the displacements from the inverse method
            writeDisplacementFF(write_options.write_path,[u_inv(1:3:end)', u_inv(2:3:end)', u_inv(3:3:end)'],'inv');
            
            %Calculate the stresses
            runFFresultFile(tools_freefem_folder_name,write_options,file_name,result_folders,code_names,mech_props,'inv');              
            
            %Gerate all the necessary mechanical variables
            write_options.write_path = result_folders.freefem_results_folder;
            mechanical_variables = generateMechanicalVariables(write_options.write_path,elem_gel,elem_cell,u_fwd,u_inv,t_fwd,t_inv,coordinates,faces,mech_props.solver);        

    end
    
    %Generate visualization files
    d.Message = sprintf(['Generating visualization files of ' file_name(1:end-4) '...']);
    
    %Generate the ouptut data and visualization files
    generateResultFiles(result_folders,write_options.write_filename,coordinates,u_fwd,u_inv,faces,mechanical_variables,visualization_tool,elem_gel,elem_cell);
    
    %Copy the template for Paraview if the selected visualization tool was
    %Paraview
    if or(strcmp(visualization_tool,'paraview'),strcmp(visualization_tool,'both'))
        copyfile([tools_paraview_folder_name filesep 'template.pvsm'],result_folders.paraview_results_folder);
    end
    disp(['Finished computing tractions for file ', file_name(1:end-4), '!!']);
    
    % Store the meshing parameters
    writeParameterFile(save_path,mesh_options,size(elem_gel,1),size(elem_cell,1),crop_rect,mech_props);    
    
    save([save_path filesep file_names.trac_param],'mech_props','mesh_options');
    
    clear u_fwd u_inv t_fwd t_inv coordinates faces;
end

% Close the dialog box
close(d);
delete(f);
