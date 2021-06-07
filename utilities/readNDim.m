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

function imOut = readNDim(fileName,info)
% This function reads multipage tiff images


mtdt = extractTifMetadata(fileName);


warning('off','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning')
warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning')


if isempty(info.t)
    mtdt.tNum = info.tNum;
    if mtdt.tNum==1
        mtdt.t = 1;
    end
else
    mtdt.tNum = 1 ; % number of images in time to be read
    mtdt.t = info.t;
end
if isempty(info.z)
    mtdt.zNum = info.zNum; 
    if mtdt.zNum==1 
        mtdt.z = 1;
    end
else
    mtdt.zNum = 1 ; % number of images in z to be read
    mtdt.z = info.z;
end
switch mtdt.bitDepth
    case 8
        dataType = 'uint8';
    case 12
        dataType = 'uint16';
    case 16
        dataType = 'uint16';
    case 24
        dataType = 'uint32';
end


if  isfield(info, 'hyperStack')
    switch  info.hyperStack
        case 'ZT'
            imOut=squeeze(zeros(mtdt.rowNum,mtdt.colNum,mtdt.zNum,mtdt.tNum,dataType));
        case 'TZ'
            imOut=squeeze(zeros(mtdt.rowNum,mtdt.colNum,mtdt.tNum,mtdt.zNum,dataType));
    end
else
    imOut=squeeze(zeros(mtdt.rowNum,mtdt.colNum,mtdt.zNum,mtdt.tNum,dataType));
end





TifLink = Tiff(fileName, 'r');

if (mtdt.tNum==1)
    tIndx =mtdt.t;
else
    tIndx = 1:mtdt.tNum;
end
if (mtdt.zNum==1)
    zIndx =mtdt.z;
else
    zIndx = 1:mtdt.zNum;
end
if (info.zNum>1) && (info.tNum>1)
    switch info.hyperStack
        case 'ZT'
            for tt=tIndx
                for zz=zIndx
                    imIndx = zz + ((tt-1)*mtdt.zNum);
                    TifLink.setDirectory(imIndx);
                    if (length(tIndx)==1)&&(length(zIndx)==1)
                        imOut(:,:)=TifLink.read();
                    elseif (length(tIndx)==1)&&(length(zIndx)>1)
                        imOut(:,:,zz)=TifLink.read();
                    elseif (length(tIndx)>1)&&(length(zIndx)==1)
                        imOut(:,:,tt)=TifLink.read();
                    elseif (length(tIndx)>1)&&(length(zIndx)>1)
                        imOut(:,:,zz,tt)=TifLink.read();
                    end
                end
            end
        case 'TZ'
            for zz=zIndx
                for tt=tIndx
                    imIndx = tt + ((zz-1)*mtdt.tNum);
                    TifLink.setDirectory(imIndx);
                    if (length(tIndx)==1)&&(length(zIndx)==1)
                        imOut(:,:)=TifLink.read();
                    elseif (length(tIndx)==1)&&(length(zIndx)>1)
                        imOut(:,:,zz)=TifLink.read();
                    elseif (length(tIndx)>1)&&(length(zIndx)==1)
                        imOut(:,:,tt)=TifLink.read();
                    elseif (length(tIndx)>1)&&(length(zIndx)>1)
                        imOut(:,:,tt,zz)=TifLink.read();
                    end
                end
            end
    end
elseif (info.zNum>1) && (info.tNum==1)
    for zz=zIndx
        imIndx = zz;
        TifLink.setDirectory(imIndx);
        if length(zIndx)==1
            imOut(:,:)=TifLink.read();
        else
            imOut(:,:,zz)=TifLink.read();
        end
    end
elseif (info.zNum==1) && (info.tNum>1)
    for tt=tIndx
        imIndx = tt ;
        TifLink.setDirectory(imIndx);
        if length(tIndx)==1
            imOut(:,:)=TifLink.read();
        else
            imOut(:,:,tt)=TifLink.read();
        end
    end
elseif (info.zNum==1) && (info.tNum==1)
    TifLink.setDirectory(1);
    imOut(:,:) = TifLink.read();
end    

    
    
 
TifLink.close();

warning('on','MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning')
warning('on','MATLAB:imagesci:tiffmexutils:libtiffWarning')


