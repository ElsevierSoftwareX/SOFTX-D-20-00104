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
function [tt_k,uu,elapsed_time_inv] = calculateNonLinInverseDisplacementField(result_folders,write_options,mech_props,elem_gel,elem_cell,coordinates,faces,solver_name,u_fwd)

%TEST!!! 
% u_fwd = u_fwd/5;

% figure
% plotmesh(coordinates,faces,[elem_cell;elem_gel]);
% daspect([1 1 1]);

%Save the file name (because we're going to modify it in each iteration)
filename = write_options.write_filename;


%Initialization
k_non = 0;
uu_trial = zeros(1,3*size(coordinates,1)); 
TOL = 1e-2;
eps = inf;
tic
while eps > TOL
    
    %Get the displacements from the previous iteration
    uu = uu_trial;
    
%     %TEST
%     if k_non == 1
%         uu = u_fwd;
%     end
    
    %Restore the original filename (for the first iteration this doesn't do
    %anything, it's for the following one, as the name changes according to
    %the iteration)
    write_options.write_filename = filename;
    
    %Write the Abaqus script    
    new_coords = coordinates-[uu(1:3:end)', uu(2:3:end)', uu(3:3:end)'];
    new_elem = [elem_gel;elem_cell];
    vol=elemvolume (new_coords(:,1:3),new_elem(:,1:4));
    q = meshquality(new_coords,new_elem);
    writeAbaqusControlStMtxScript(new_coords,elem_gel,elem_cell,write_options,mech_props,k_non,uu);
  
    disp(['Iteration ' num2str(k_non) '... Obtaining the stiffness matrix from Abaqus...']);
    %Change the name and run the script
    write_options.write_filename = [filename '_' num2str(k_non)];
    runAbaqusFile(write_options,mech_props,'stiffness_matrix');
    
    %Load the stiffness matrix
    A = getStiffnessMatrix(result_folders,write_options,mech_props,[elem_gel;elem_cell],coordinates);
    
    %Load the forces
    tt_aux = load([result_folders.intermediate_results_folder filesep 'reac_for.dat']);
    tt_k = reshape(tt_aux(:,2:end)',3*size(coordinates,1),1);
    
    %Build the matrix system of equations
    ind = unique([unique(elem_cell(:,1:4));unique(faces(:,1:3))]);
    ind_a = [3*ind-2;3*ind-1;3*ind]; clear ind;
    ind_b = setdiff(1:length(u_fwd),ind_a);
    ua_s = u_fwd(ind_a)';
    ub_s = u_fwd(ind_b)';
    b_eval = -tt_k(ind_b) + A(ind_b,ind_a)*uu(ind_a)'+A(ind_b,ind_b)*uu(ind_b)';
    ia = length(ind_a);
    ib = length(ind_b);
    
    % save(strcat('b_eval_',num2str(k_non)),'b_eval');
    save([result_folders.intermediate_results_folder filesep filename '_tt_k_',num2str(k_non)],'tt_k');
    
    switch solver_name
        case 'gmres'
            x = gmres(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; b_eval],[],1e-3,5500);
        case 'suitesparse'
            x = spqr_solve([speye(ia) speye(ia,ib)-speye(ia,ib) A(ind_a,ind_b);...
                speye(ib,ia)-speye(ib,ia) speye(ib) A(ind_b,ind_b);...
                A(ind_b,ind_a) A(ind_b,ind_b) speye(ib,ib)-speye(ib,ib)],...
                [ua_s; ub_s; b_eval]);
    end
    
    uu_trial(ind_a) = x(1:length(ind_a));
    uu_trial(ind_b) = x(length(ind_a)+1:length(ind_a)+length(ind_b));
    
    eps = (norm(uu_trial-uu))/norm(uu_trial)
    
    k_non = k_non+1        

%     save(strcat('uu_trial_',num2str(k_non)),'uu_trial');
    
end
elapsed_time_inv = toc;
