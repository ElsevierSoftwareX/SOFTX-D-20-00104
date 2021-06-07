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

function data = diffSurfCalc(diffcellmask,resolution)


% Prepare the coordinates and tringulation appropriately
pointCount = 0;
data.vertices = [];
data.faces = [];
data.colorIndx = [];
for ii=1:2
    if ii==1 
        mask = diffcellmask>0; % protrusions
    elseif ii==2
        mask = diffcellmask<0; % retractions
    end 
    % Smooth slightly the mask
    blurredMask = im2mat(gaussf(mask,1)); clear mask
    % triangulation of the surface
    th = 0.1;
    [fv,~] = getSurfNormals(blurredMask,th,resolution);
    pointNum = size(fv.vertices,1);
    if pointNum>0
        data.colorIndx = [data.colorIndx; repmat(ii-1,pointNum,1)];
        data.vertices = [data.vertices;fv.vertices];
        data.faces = [data.faces; fv.faces + pointCount];
        pointCount = pointCount + pointNum;
    end
    clear fv blurredMask pointNum
end


