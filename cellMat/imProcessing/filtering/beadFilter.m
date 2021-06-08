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

function out = beadFilter(im,var1,var2,weight,smoothVar,lb,ub)
% This function applies a filtering step to the bead images. In some cases,
% the filtering step can remove the background light produced by widefield
% (non-confocal) microscopes
% The filtering combines the Difference-of-Gaussian operator with some
% contrast stretching

%% No DIPImage
imDim =  length(size(im));

% Pad
padVal = 5;
im = padarray(double(im),padVal*ones(1,imDim),'replicate');

% Filter
im1 = imgaussfilt3(im,var1);
im2 = weight*imgaussfilt3(im,var2);
im = im1-im2;
clear im1 im2;
im = imadjustn(mat2gray(im),stretchlim(mat2gray(im(:)),[lb ub]/100),[0,1])*255;
if smoothVar>0
    im = imgaussfilt3(im,smoothVar);
end
im = imadjustn(mat2gray(im),[0 1],[0,1])*255;
% dip_image(im)

% Remove pad
switch imDim
    case 2
        out = im(padVal+1:end-padVal,padVal+1:end-padVal);
    case 3
        out = im(padVal+1:end-padVal,padVal+1:end-padVal,padVal+1:end-padVal);
end
clearvars -except out;

out = uint8(round(out));

%% Using DIPImage functions:
% imDim =  length(size(im));
% 
% % Pad
% padVal = 5;
% im = dip_image(padarray(double(im),padVal*ones(1,imDim),'replicate'));
% % Filter
% out = stretch(gaussf(stretch(gaussf(im,var1)- weight*gaussf(im,var2),lb,ub,0,255),smoothVar));
% % Remove pad
% switch imDim
%     case 2
%         out = out(padVal:end-padVal,padVal:end-padVal);
%     case 3
%         out = out(padVal:end-padVal,padVal:end-padVal,padVal:end-padVal);
% end
% 
% out = im2mat(round(out),'uint8');

%% PROBATINA!
% %Denoising
% out = medfilt3(out);
% %Contrast enhancement
% out = imadjustn(out);



