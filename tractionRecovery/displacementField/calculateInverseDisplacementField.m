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
function [u_inv,elapsed_time_inv] = calculateInverseDisplacementField(elem,fc1,u_fwd,A,solver_name)

%Build the matrix system of equations
ind = unique([unique(elem(:,1:4));unique(fc1(:,1:3))]);
ind_a = [3*ind-2;3*ind-1;3*ind]; clear ind;
ind_b = setdiff(1:length(u_fwd),ind_a);
ua_s = u_fwd(ind_a)';
ub_s = u_fwd(ind_b)';
ia = length(ind_a);
ib = length(ind_b);

tic
switch solver_name    
    case 'gmres'        
        %Solve using Matlab function
        x = gmres(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; zeros(length(ind_b),1)],[],1e-3,30000);
        elapsed_time_inv = toc;
    case 'suitesparse'
        %Solve using SuiteSparse
%         bb = [speye(ia) speye(ia,ib)-speye(ia,ib) A(ind_a,ind_b); speye(ib,ia)-speye(ib,ia) speye(ib) A(ind_b,ind_b); A(ind_b,ind_a) A(ind_b,ind_b) speye(ib,ib)-speye(ib,ib)];
%         AA = [ua_s; ub_s; zeros(length(ind_b),1)];
        tic
        x = spqr_solve([speye(ia) speye(ia,ib)-speye(ia,ib) A(ind_a,ind_b); speye(ib,ia)-speye(ib,ia) speye(ib) A(ind_b,ind_b); A(ind_b,ind_a) A(ind_b,ind_b) speye(ib,ib)-speye(ib,ib)],...
            [ua_s; ub_s; zeros(length(ind_b),1)]);
        elapsed_time_inv = toc;
end

%Build the displacement field
u_inv(ind_a) = x(1:length(ind_a));
u_inv(ind_b) = x(length(ind_a)+1:length(ind_a)+length(ind_b));
