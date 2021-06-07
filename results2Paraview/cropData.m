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

function data = cropData(data,roi)

fieldSize = size(data);

if isempty(roi.rowRange)
    rowRange = [1,fieldSize(1)];
else
    rowRange = roi.rowRange;
end
if isempty(roi.colRange)
    colRange = [1,fieldSize(2)];
else
    colRange = roi.colRange;
end
if isempty(roi.depthRange)
    depthRange = [1,fieldSize(3)];
else
    depthRange = roi.depthRange;
end
data = squeeze(data(rowRange(1):rowRange(2),colRange(1):colRange(2),depthRange(1):depthRange(2)));

end