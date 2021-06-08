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


function dispCalc(timePoints,pathParam,dispParam)
% This function calculates the displacements in the timelapse using
% FFD-based image registration. Different strategies can be used.
% Then the calculated displacementes are corrected against different
% artifacts



% Calcualte the displacements
timeSeriesDispCalc(timePoints.pending,pathParam,dispParam);


% Shift correction (translational motion) <-- if required
 if dispParam.shiftCorr
     timeSeriesDispShiftCorr(timePoints.all,pathParam,dispParam.strategy);
 end


% Correction of artifacts
if dispParam.artifactCorr
    timeSeriesDispArtifactCorr(timePoints.all,pathParam) 
end
