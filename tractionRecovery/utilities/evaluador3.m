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

function f = evaluador(x,A,ind_a,ind_b)


ib = length(ind_b);
ia = length(ind_a);


a1 = x(1:ia)+A(ind_a,ind_b)*x(ia+ib+1:end);

a2 = x(ia+1:ia+ib)+A(ind_b,ind_b)*x(ia+ib+1:end);

a3 = A(ind_b,ind_a)*x(1:ia)+A(ind_b,ind_b)*x(ia+1:ia+ib);




f = [a1; a2; a3];

