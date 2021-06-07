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
function [C10,D1] = getNeoHookeanParameters(E,v)

%Relations: https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/calc_elastic_constants.cfm#calc

%Shear modulus (mu, G)
mu = E / (2*(1+v));

%Bulk modulus (k, K)
k = E / (3*(1-2*v));

%Calculate the temperature dependent material parameters
C10 = mu/2;
D1 = 2/k;
