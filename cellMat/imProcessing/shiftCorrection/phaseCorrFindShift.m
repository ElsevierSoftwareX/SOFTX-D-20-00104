
%    This file is part of cellMat.
%     Copyright (C) 2016 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
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


function [shft,xcPeak] = phaseCorrFindShift(in1,in2)

dim = length(size(in1));

ft1 = fftn(in1);
ft2 = fftn(in2);

% Phase Correlation 
XC = ft1.*conj(ft2); 
XC = XC./abs(XC); 
xc = real(ifftshift(ifftn(XC)));
fieldSize = size(xc);
clear XC ft2

% Phase autocorrelation
AC = ft1.*conj(ft1);
AC = AC./abs(AC);
ac = real(ifftshift(ifftn(AC)));
clear AC ft1


% Shift calc

switch dim
    case 2
        [xcPeak,xcIndx] = max(xc(:));
        [xc_peakRow,xc_peakCol] = ind2sub(fieldSize,xcIndx);
        clear xc xcIndx
        [~,acIndx] = max(ac(:));
        [ac_peakRow,ac_peakCol] = ind2sub(fieldSize,acIndx);
        clear ac acIndx
        shft = [xc_peakRow,xc_peakCol] - [ac_peakRow,ac_peakCol];
        shft = [shft(2),shft(1)]; % arranged as [x,y]
    case 3
        [xcPeak,xcIndx] = max(xc(:));
        [xc_peakRow,xc_peakCol,xc_peakDepth] = ind2sub(fieldSize,xcIndx);
        clear xc xcIndx
        [~,acIndx] = max(ac(:));
        [ac_peakRow,ac_peakCol,ac_peakDepth] = ind2sub(fieldSize,acIndx);
        clear ac acIndx
        shft = [xc_peakRow,xc_peakCol,xc_peakDepth] - [ac_peakRow,ac_peakCol,ac_peakDepth];
        shft = [shft(2),shft(1),shft(3)]; % arranged as [x,y,z]
end





    