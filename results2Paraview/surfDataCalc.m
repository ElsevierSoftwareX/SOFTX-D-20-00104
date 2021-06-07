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

function data = surfDataCalc(mask,vecField,resolution)

% Smooth slightly the mask
blurredMask = im2mat(gaussf(mask,1));

% Triangulation of the cell surface and corresponding normals
th =0.1;
filter_flag = true;
[data,surfNormVec] = getSurfNormals(blurredMask,th,resolution,filter_flag);

% vector field at the surface
surfVecField = vectorProjSurf(vecField,data.vertices,resolution);

% vecotor field magnitude at the surface
data.vecMag = sqrt(surfVecField.X(:).^2 + surfVecField.Y(:).^2 + surfVecField.Z(:).^2);

% vector field normal to the surface
surfNormVec = surfNormVec./repmat(sqrt(sum(surfNormVec.^2,2)),[1,3]); % Make sure that the vector normal to the surface (surfNorm) es normalized
data.vecNormComp = surfNormVec(:,1).*surfVecField.X(:) + surfNormVec(:,2).*surfVecField.Y(:) + surfNormVec(:,3).*surfVecField.Z(:);

% vector field parallel to the surface
surfShearVec  = [surfVecField.X(:),surfVecField.Y(:),surfVecField.Z(:)] - (repmat(data.vecNormComp,[1,3]).*surfNormVec);
data.vecShearComp = sqrt(sum(surfShearVec.^2,2));

% Normal and parallel components percentage
tot = abs(data.vecNormComp) + data.vecShearComp;
data.vecNormPercent = 100*(data.vecNormComp./tot);
data.vecNormPercent(isnan(data.vecNormPercent)) = 0;
data.vecShearPercent = 100*(data.vecShearComp./tot);
data.vecShearPercent(isnan(data.vecShearPercent)) = 0;
