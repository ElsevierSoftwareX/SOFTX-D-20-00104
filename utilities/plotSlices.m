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

function plotSlices(handle_axesOriginal,handle_axesFiltered,channel,slice)
im_noFilt = getappdata(0,['im_',channel,'_noFilt']);
im_filt = getappdata(0,['im_',channel,'_filt']);
imshow(im_noFilt(:,:,round(slice)),[],'Parent',handle_axesOriginal);
if isempty(im_filt)
    %Plot nothing
    imshow(im_filt,[],'Parent',handle_axesFiltered);
else 
    %Plot the selected slice
    if strcmp(channel,'cellSeg')
        tmp = labeloverlay(im_noFilt(:,:,round(slice)),im_filt(:,:,round(slice)),'Transparency',0.8);
        imshow(tmp,[],'Parent',handle_axesFiltered);        
    else
        imshow(im_filt(:,:,round(slice)),[],'Parent',handle_axesFiltered);
    end
end