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

function [outputVec,normalComp,shearComp,shearVec] = tensorProjSurf(input,surfNorm,surfVert,mask,extrapFun,resolution)

% It is assumed that input is a square tensor of size 3x3

fieldSize = size(input{1,1});
dim = 3;
surfVertNum = size(surfVert,1);


% Make sure that the vector normal to the surface (surfNorm) es normalized
surfNorm = surfNorm./repmat(sqrt(sum(surfNorm.^2,2)),[1,dim]);

% Extrapolate the input tensor to the vertices of the 
[rCoord,cCoord,dCoord] = find(mask);
padVal = 5;
calcMask = false(fieldSize);
calcMask(max(1,min(rCoord)-padVal):min(fieldSize(1),max(rCoord)+padVal), ...
         max(1,min(cCoord)-padVal):min(fieldSize(2),max(cCoord)+padVal), ...
         max(1,min(dCoord)-padVal):min(fieldSize(3),max(dCoord)+padVal)) = true;
calcMask(mask(:)>0) = false;

inputSurf = cell([dim,1]);
[inputSurf{:}]=deal(zeros([surfVertNum,1]));
for ii=1:dim
    for jj=1:dim
        extrapFun.Values = input{ii,jj}(calcMask(:));
        inputSurf{ii,jj} = extrapFun(surfVert(:,2),surfVert(:,1),surfVert(:,3));
    end
end


% Main stuff
outputVec = zeros([surfVertNum,dim]);

outputVec(:,1) = surfNorm(:,1).*inputSurf{1,1} + surfNorm(:,2).*inputSurf{1,2} + surfNorm(:,3).*inputSurf{1,3};
outputVec(:,2) = surfNorm(:,1).*inputSurf{2,1} + surfNorm(:,2).*inputSurf{2,2} + surfNorm(:,3).*inputSurf{2,3};
outputVec(:,3) = surfNorm(:,1).*inputSurf{3,1} + surfNorm(:,2).*inputSurf{3,2} + surfNorm(:,3).*inputSurf{3,3};
normalComp = surfNorm(:,1).*outputVec(:,1) + surfNorm(:,2).*outputVec(:,2) + surfNorm(:,3).*outputVec(:,3);
shearVec  = outputVec - (repmat(normalComp,[1,dim]).*surfNorm);
shearComp = sqrt(sum(shearVec.^2,2));
shearVec = shearVec./repmat(shearComp,[1,dim]); % normalization
% Alternative way of getting just the magnitude of the sehar component:
    % outputNorm2 = (outputVec(:,1).^2 + outputVec(:,2).^2 + outputVec(:,3).^2); % squared norm
    % shearComp = sqrt(outputNorm2 - normalComp.^2);




