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

function TF = tractionforces(coordinates,fc1,stress,elem0)

%Compute the direction vector for a segment of the element face
s1_aux = coordinates(fc1(:,2),:)-coordinates(fc1(:,1),:);

%Compute the direction vector for another segment of the element face
s2_aux = coordinates(fc1(:,3),:)-coordinates(fc1(:,2),:);

%Compute the perpendicular vectors to the element face
sur = cross(s1_aux,s2_aux,2);

%Remove the useless column of the element matrix and the faces matrix
elem0_cut = elem0(:,1:4);
fc1_cut = fc1(:,1:3);

%Sort the lists of nodes in order
elem0_sort = sort(elem0_cut,2);
fc1_sort = sort(fc1_cut,2);


%Look for the faces that the ECM and the cell have in common:
[~, I1] = ismember(fc1_sort,elem0_sort(:,[1 2 3]),'rows');
[~, I2] = ismember(fc1_sort,elem0_sort(:,[1 3 4]),'rows');
[~, I3] = ismember(fc1_sort,elem0_sort(:,[1 2 4]),'rows');
[~, I4] = ismember(fc1_sort,elem0_sort(:,[2 3 4]),'rows');

%ismember returns zero when the face is not found, so let's check where the
%zeros are
[i_q,j_q] = find([I1 I2 I3 I4]~=0);
aux_q = sortrows([i_q j_q]);
clear i_q j_q;

%Gather all the indices in one vector by summing them 
%(so that there are no zeros left)
IA = sum([I1 I2 I3 I4],2);

%Store the columns that contain zeros
j_qq = aux_q(:,2);

p4_aux = [4 2 3 1];
in_q = sub2ind([size(elem0_cut,1),size(elem0_cut,2)],IA,p4_aux(j_qq)');

%Get all the indices of the common nodes between cell and ECM
p4 = elem0_sort(in_q);

%Correct the sign of the normal vectors
v4=coordinates(p4,:)-coordinates(fc1(:,1),:);
sf2=sign(sum(v4.*sur,2));

suri_x = sf2.*sur(:,1).*ones(size(fc1,1),3)/3;
suri_y = sf2.*sur(:,2).*ones(size(fc1,1),3)/3;
suri_z = sf2.*sur(:,3).*ones(size(fc1,1),3)/3;


Igg = fc1(:,[1,2,3]);
surii_x = sparse ( Igg , ones(size(fc1,1),3) , suri_x, size(coordinates,1) , 1 ) ; 
surii_y = sparse ( Igg , ones(size(fc1,1),3) , suri_y, size(coordinates,1) , 1 ) ; 
surii_z = sparse ( Igg , ones(size(fc1,1),3) , suri_z, size(coordinates,1) , 1 ) ; 

surii_n = sqrt(surii_x.*surii_x + surii_y.*surii_y + surii_z.*surii_z);

nor_sx(1:size(coordinates,1),1) = 0;
nor_sy(1:size(coordinates,1),1) = 0;
nor_sz(1:size(coordinates,1),1) = 0;

nor_sx(unique(fc1(:,1:3))) = surii_x(unique(fc1(:,1:3)))./surii_n(unique(fc1(:,1:3)));
nor_sy(unique(fc1(:,1:3))) = surii_y(unique(fc1(:,1:3)))./surii_n(unique(fc1(:,1:3)));
nor_sz(unique(fc1(:,1:3))) = surii_z(unique(fc1(:,1:3)))./surii_n(unique(fc1(:,1:3)));

numnode = 1:size(coordinates,1);


RF1 = stress(:,1).*nor_sx + stress(:,4).*nor_sy + stress(:,5).*nor_sz;
RF2 = stress(:,4).*nor_sx + stress(:,2).*nor_sy + stress(:,6).*nor_sz;
RF3 = stress(:,5).*nor_sx + stress(:,6).*nor_sy + stress(:,3).*nor_sz;

TF(3*numnode-2) = -RF1;
TF(3*numnode-1) = -RF2;
TF(3*numnode) = -RF3;

% TF(3*numnode-2) = nor_sx;
% TF(3*numnode-1) = nor_sy;
% TF(3*numnode) = nor_sz;