%    This file is part of TFMLAB.
%     Copyright (C) 2016-2020 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
%     Copyright (C) 2019-2020 Jorge Barrasa-Fano
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


function [shft,ax_ang] = shiftCalc(imA,imB,options,file_name)
% This function calculates the shift (rigid translation) betwen two images


if not(strcmp(options.method,'Translation phaseCorr (fast)')) || not(strcmp(options.method,'Translation Elastix (medium)')) || not(strcmp(options.method,'Euler Elastix (slow)'))
    % sigma for the gaussina smooting when computing the shift
    shiftSmoothSigma = 3; 
    % Change the boundary condtions to padd with zeros the locations where the
    % image won't exist due to the shift correction
    oldBoundaryOption = dip_getboundary(1);
    dip_setboundary('add_zeros');
end

ax_ang = [];
switch options.method
    case 'Euler Elastix (slow)'
        options.elastix.txType = 'EulerTransform';
        [shft,ax_ang] = elastixFindShift(imA,imB,options.elastix,file_name);
    case 'Translation phaseCorr (fast)'
        shft = phaseCorrFindShift(imA,imB);
    case 'Translation Elastix (medium)'
        options.elastix.txType = 'TranslationTransform';
        [shft,ax_ang] = elastixFindShift(imA,imB,options.elastix,file_name);
    otherwise
        shft = findshift(dip_image(imA),dip_image(imB),options.method,shiftSmoothSigma)';
end




if not(strcmp(options.method,'elastix'))||not(strcmp(options.method,'phaseCorr'))
    % Change the boundary conditions back to its previous value
    dip_setboundary(oldBoundaryOption);
end