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
function [node_all,el1_10]= getC3D10elements(node,el1)

%Coordinates of the new nodes 5-10
nodes5 = (node(el1(:,1),:)+ node(el1(:,2),:))/2; 
nodes6 = (node(el1(:,2),:)+ node(el1(:,3),:))/2;
nodes7 = (node(el1(:,1),:)+ node(el1(:,3),:))/2;
nodes8 = (node(el1(:,1),:)+ node(el1(:,4),:))/2; 
nodes9 = (node(el1(:,2),:)+ node(el1(:,4),:))/2; 
nodes10 = (node(el1(:,3),:)+ node(el1(:,4),:))/2;
new_nodes = [nodes5;nodes6;nodes7;nodes8;... 
    nodes9;nodes10]; 
%Clear all the individual nodes
clear nodes*, 

% C = A(ia,:) and A = C(ic,:).
[C,~,ic] = unique(new_nodes,'rows','stable');

%tmp(ic) gives me a label for each node
tmp = 1:length(C);
%I add a new column that gives me the node identifier.
new_nodes = [new_nodes tmp(ic)'];
clear tmp ic;

%Join all the nodes
node_all = [node;C];
clear C;

%Update the labels:
new_nodes(:,4) = new_nodes(:,4)+length(node);

%Calculate the number of elements
nEls = length(el1);

%Create a new variable for the 10 node elements
el1_10 = zeros(nEls,11);

%Store the ones that I had before
el1_10(:,1:4) = el1(:,1:4);
el1_10(:,11) = el1(:,5);

%Store the new ones
el1_10(:,5)  = new_nodes(1       :1*nEls,4);
el1_10(:,6)  = new_nodes(1*nEls+1:2*nEls,4);
el1_10(:,7)  = new_nodes(2*nEls+1:3*nEls,4);
el1_10(:,8)  = new_nodes(3*nEls+1:4*nEls,4);
el1_10(:,9)  = new_nodes(4*nEls+1:5*nEls,4);
el1_10(:,10) = new_nodes(5*nEls+1:6*nEls,4);

% 
% figure
% plotmesh(node_all,faces)
% hold on;
% scatter3(node_all(:,1),node_all(:,2),node_all(:,3),8,'g','filled');






