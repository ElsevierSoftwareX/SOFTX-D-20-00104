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


function [fv,n] = getSurfNormals(varargin)

im  = varargin{1};
th  = varargin{2};
resolution  = varargin{3};

fieldSize = size(im);
[x,y,z] = meshgrid(1:fieldSize(2),1:fieldSize(1),1:fieldSize(3));
x = x*resolution.xy;
y = y*resolution.xy;
z = z*resolution.z;

fv = isosurface(x,y,z,im,th);
if nargin > 3 && varargin{4}
    filter_options.mode = 1;
    filter_options.itt = 7;
    filter_options.lambda = 1;
    filter_options.sigma = 1;
    try
        fv=smoothpatch(fv,filter_options.mode,filter_options.itt,filter_options.lambda,filter_options.sigma);
    catch exception
        fv=smoothpatch(fv,filter_options.mode,filter_options.itt,filter_options.lambda,filter_options.sigma);
    end    
end
n = isonormals(x,y,z,im,fv.vertices); % sometimes it provides a weird error

