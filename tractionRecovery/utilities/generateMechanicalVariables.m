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
function mechanical_variables = generateMechanicalVariables(results_path,elem0,elem,u_fwd,u_inv,t_fwd,t_inv,coordinates,faces,solver)
%Group all the elements in one matrix
elements = [elem0;elem];

%Compute volume forces
BDF_i = volumeforces(coordinates,elements,t_inv);
BDF_f = volumeforces(coordinates,elements,t_fwd);

%Create a container for the stresses and strains
stress_f(1:size(coordinates,1),1:6) = 0;
stress_i(1:size(coordinates,1),1:6) = 0;
stressp_f(1:size(coordinates,1),1:3) = 0;
stressp_i(1:size(coordinates,1),1:3) = 0;
strainp_f(1:size(coordinates,1),1:3) = 0;
strainp_i(1:size(coordinates,1),1:3) = 0;

% Load the data from Abaqus
switch solver
    case 'Abaqus'
        load([results_path filesep 'stresses.dat']);
        load([results_path filesep 'stressesp.dat']);
        load([results_path filesep 'strainsp.dat']);        
        
        stress_f(unique(elem0(:,1:4)),:) = stresses(1:length(unique(elem0(:,1:4))),2:end);
        stress_i(unique(elem0(:,1:4)),:) = stresses(length(unique(elem0(:,1:4)))+1:end,2:end);
        
        stressp_f(unique(elem0(:,1:4)),:) = stressesp(1:length(unique(elem0(:,1:4))),2:end);
        stressp_i(unique(elem0(:,1:4)),:) = stressesp(length(unique(elem0(:,1:4)))+1:end,2:end);
        
        strainp_f(unique(elem0(:,1:4)),:) = strainsp(1:length(unique(elem0(:,1:4))),2:end);
        strainp_i(unique(elem0(:,1:4)),:) = strainsp(length(unique(elem0(:,1:4)))+1:end,2:end);
        
    case 'FreeFEM'
        load([results_path filesep 'stresses_fwd.txt']);
        load([results_path filesep 'stressesp_fwd.txt']);
        load([results_path filesep 'strainsp_fwd.txt']);
        load([results_path filesep 'stresses_inv.txt']);
        load([results_path filesep 'stressesp_inv.txt']);
        load([results_path filesep 'strainsp_inv.txt']);
        
%         indx_nodes_gel = unique(elem0(:,1:4));
        selected_indx = ismember(stresses_fwd(:,1),unique(elem0(:,1:4)));
        
        stress_f(selected_indx,:) = stresses_fwd(selected_indx,3:end);
        stress_i(selected_indx,:) = stresses_inv(selected_indx,3:end);
        stressp_f(selected_indx,:) = stressesp_fwd(selected_indx,3:end);
        stressp_i(selected_indx,:) = stressesp_inv(selected_indx,3:end);
        strainp_f(selected_indx,:) = strainsp_fwd(selected_indx,3:end);
        strainp_i(selected_indx,:) = strainsp_inv(selected_indx,3:end);
end

TF_i = tractionforces(coordinates,faces,stress_i,elem0);
TF_f = tractionforces(coordinates,faces,stress_f,elem0);

%Build a result container
mechanical_variables.uf = u_fwd';
mechanical_variables.ui = u_inv';
mechanical_variables.bf =BDF_f';
mechanical_variables.bi =BDF_i';
mechanical_variables.tf =TF_f';
mechanical_variables.ti =TF_i';
mechanical_variables.ssf =stressp_f;
mechanical_variables.ssi =stressp_i;
mechanical_variables.srf =strainp_f;
mechanical_variables.sri =strainp_i;