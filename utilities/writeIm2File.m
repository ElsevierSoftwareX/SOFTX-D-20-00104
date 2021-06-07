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

function writeIm2File(im,path,fileName,imFormat)
% This function writes a 2D/3D image to a file

dim = length(size(im));


switch imFormat
    
    case 'tiff'
        % Default resolution value: 72 px/inch (for 2D)
        switch dim
            case 2
                imRes = 72; % 72 px/inch
                imwrite(im,[path filesep fileName '.' imFormat],imFormat,'Resolution',imRes);
            case 3
                if not(exist([path filesep fileName '.' imFormat],'file'))
                    imwrite(squeeze(im(:,:,1)),[path filesep fileName '.' imFormat],imFormat);
                    for ii=2:size(im,3)
                        imwrite(squeeze(im(:,:,ii)),[path filesep fileName '.' imFormat],imFormat,'WriteMode','append');
                    end
                else
                    for ii=1:size(im,3)
                        imwrite(squeeze(im(:,:,ii)),[path filesep fileName '.' imFormat],imFormat,'WriteMode','append');
                    end
                end
                
        end
    case 'nii' % TO BE CHECKED --> RESOLUTION UNIT AND VALUE
        niiIm = make_nii(im); clear im
        save_nii(niiIm, [path filesep fileName '.' imFormat], false);
end

