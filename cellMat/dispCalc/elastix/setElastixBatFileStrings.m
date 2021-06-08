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

function [rmv_line,str_in,str_out] = setElastixBatFileStrings(elastixLibFolder,im_folder,param_folder,fixed_im_file_name,moving_im_file_name,im_format,fixed_mask_name,moving_mask_name,dispParam)

count = 1;
count_rmv =1; 
rmv_line{count_rmv} = [];


str_in{count} = 'elastixLibPath'; str_out{count} = elastixLibFolder; count=count+1;
str_in{count} = 'imReadPath'; str_out{count} = im_folder; count=count+1;
str_in{count} = 'parametersReadPath'; str_out{count} = param_folder; count=count+1;
str_in{count} = 'fixedImName'; str_out{count} = [fixed_im_file_name '.' im_format];  count=count+1;
str_in{count} = 'movingImName'; str_out{count} = [moving_im_file_name '.' im_format]; count=count+1;
str_in{count} = 'parameterFileName'; str_out{count} = 'elastixDispCalcParam.txt'; count=count+1;
if dispParam.mask.flag
    str_in{count} = 'fixedMaskName'; str_out{count} = [fixed_mask_name '.' im_format]; count=count+1;
    str_in{count} = 'movingMaskName'; str_out{count} = [moving_mask_name '.' im_format]; count=count+1;
else
    rmv_line{count_rmv} = 'SET fixedMask=%imPath%\fixedMaskName'; count_rmv = count_rmv+1;
    rmv_line{count_rmv} = 'SET movingMask=%imPath%\movingMaskName'; count_rmv = count_rmv+1;
    rmv_line{count_rmv} = 'ECHO Fixed Mask: %fixedMask%'; count_rmv = count_rmv+1;
    rmv_line{count_rmv} = 'ECHO Moving Mask: %movingMask%'; count_rmv = count_rmv+1;
    str_in{count} = '-fMask %fixedMask% -mMask %movingMask%'; str_out{count} = '';  count=count+1;
end
if not(strcmp(dispParam.strategy,'bwdIniTf'))
    rmv_line{count_rmv} = 'SET iniTransfPath=iniTransfReadPath'; count_rmv = count_rmv+1;
    rmv_line{count_rmv} = 'SET iniTransf=%iniTransfPath%\iniTransfFileName'; count_rmv = count_rmv+1;
    str_in{count} = '-t0 %iniTransf%'; str_out{count} = ''; count=count+1;
end