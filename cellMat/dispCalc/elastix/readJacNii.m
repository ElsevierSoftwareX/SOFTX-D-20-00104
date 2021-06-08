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

function jacobMat = readJacNii(fileName1,scale,resolution)


%% Read nii data 

niiData = load_nii(fileName1);
data = squeeze(niiData.img); clear niiData

dim = length(size(data));

%% Extract the jacobian

switch dim
    case 3 % 2D
        jacobMat = cell(2,2);
    case 4 % 3D
        jacobMat = cell(3,3);
end
% [ d11 d21 d31 d12 d22 d32 d13 d23 d33] (in 3D)
% Rotate approprietely the XY-plane 
% figure
for ii = 1:numel(jacobMat)
    %     [row,col]=ind2sub([3 3],ii);
    switch dim        
        case 3
            jacobMat{ii} = flip(rot90((data(:,:,ii)),-1),1);
        case 4 % 3D
            jacobMat{ii} = flip(rot90((data(:,:,:,ii)),-1),1);
    end
    
    %     jacobMat{ii} = flip(data(:,:,:,ii),1);
    
    %     subplot(3,3,ii)
    %     imshow(flip(jacobMat{ii}(:,:,20)),[])
    %     colormap(jet(256));colorbar;
    %     title(['Component ',num2str(ii),' of Tx']);
end
%Substract the identity matrix to get the gradient tensor:
% jacobMat{1,1} = jacobMat{1,1} - 1;
% jacobMat{2,2} = jacobMat{2,2} - 1;
% if dim == 4
%     jacobMat{3,3} = jacobMat{3,3} - 1;
% end
% % Rearrange the tensor
% tmp=jacobMat;
% 
% jacobMat{1,1} = tmp{2,1};
% jacobMat{1,2} = tmp{2,2}-1;
% 
% jacobMat{2,1} = tmp{1,1}-1;
% jacobMat{2,2} = tmp{1,2};
% 
% if dim == 4 % 3D
%     jacobMat{2,3} = tmp{1,3};
%     jacobMat{1,3} = tmp{2,3};
%     jacobMat{3,3} = tmp{3,3}-1;
% end
% clear tmp
% jacobMat = jacobMat';
%% Convert the displacements to the right units
% 
%Scale the resolutions
resX = resolution.xy/scale(1);
resY = resolution.xy/scale(2);

jacobMat{1,1} = jacobMat{1,1} * (resX / resX);
jacobMat{1,2} = jacobMat{1,2} * (resX / resY);

jacobMat{2,1} = jacobMat{2,1} * (resY / resX);
jacobMat{2,2} = jacobMat{2,2} * (resY / resY);


if dim == 4
    %Scale the resolutions
    resZ = resolution.z/scale(3);
    
    jacobMat{1,3} = jacobMat{1,3} * (resX / resZ);
    
    jacobMat{2,3} = jacobMat{2,3} * (resY / resZ);
    
    jacobMat{3,1} = jacobMat{3,1} * (resZ / resX);
    jacobMat{3,2} = jacobMat{3,2} * (resZ / resY);
    jacobMat{3,3} = jacobMat{3,3} * (resZ / resZ);
end

