
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



function timeSeriesDispArtifactCorr(timePoints,pathParam)


tpNum = length(timePoints);



disp(' ... ... ...')
disp(' ... ... ...')

for tt=1:tpNum
    
    tp = timePoints(tt);
    currentFileName = [pathParam.saveFileName '_t' num2str(tp)];
    
    disp(' ... ... ...')
    disp(['Disp Artifact Corr : TIME POINT ' num2str(tp) ' ... ...'])
    
    tmp = load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'dispField');
    dispFieldArtif = tmp.dispField; clear tmp
    save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','dispFieldArtif','-v7.3')
    
    load([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'dispField');
    dispField = dFieldCorr(dispField);
    save([pathParam.savePath filesep pathParam.saveSubFolder filesep currentFileName '.mat'],'-append','dispField','-v7.3')
    clear dispField dispFieldArtif
    
end

disp('ok')

end


function dField = dFieldCorr(dField)

if  length(squeeze(size(dField.X)))==3
    s = size(dField.X);
    for ii=1:s(3)
        tmp = squeeze(dField.X(:,:,ii));
        dField.X(:,:,ii) = tmp - mean(tmp(:)); clear tmp
        tmp = squeeze(dField.Y(:,:,ii));
        dField.Y(:,:,ii) = tmp - mean(tmp(:)); clear tmp
        if isfield(dField,'Z')
            tmp = squeeze(dField.Z(:,:,ii));
            dField.Z(:,:,ii) = tmp - mean(tmp(:)); clear tmp
        end
    end
end

end

