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
function BF = volumeforces(coordinates,elements,t0);


RFX = t0(1:3:length(t0));
RFY = t0(2:3:length(t0));
RFZ = t0(3:3:length(t0));

a1_aux = coordinates(elements(:,1),:)-coordinates(elements(:,4),:);
a2_aux = coordinates(elements(:,2),:)-coordinates(elements(:,4),:);
a3_aux = coordinates(elements(:,3),:)-coordinates(elements(:,4),:);


vol = abs(sum(a1_aux.*cross(a2_aux,a3_aux,2),2))/6;




voli = vol.*ones(size(elements,1),4)/4;

Igg = elements(:,[1,2,3,4]);
volii = sparse ( Igg , ones(size(elements,1),4) , voli, size(coordinates,1) , 1 ) ;

numnode = 1:size(coordinates,1);

RFXi = RFX./volii;
RFYi = RFY./volii;
RFZi = RFZ./volii;


BF(3*numnode-2) = RFXi;
BF(3*numnode-1) = RFYi;
BF(3*numnode) = RFZi;

% BF(3*numnode-2) = volii;
% BF(3*numnode-1) = volii;
% BF(3*numnode) =  volii;

