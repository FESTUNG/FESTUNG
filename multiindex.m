% Computes all two-dimensional multi-indices involved in the representation
% of a polynomial solution of degree p.

%===============================================================================
%> @file multiindex.m
%>
%> @brief Computes all two-dimensional multi-indices involved in the representation
%>        of a polynomial solution of degree p.
%===============================================================================
%>
%> @brief  Computes all two-dimensional multi-indices involved in the representation
%>         of a polynomial solution of degree @f$p@f$.
%>
%> This can be used to determine the derivatives associated with a degree
%> of freedom in a Taylor basis representation (cf.
%> <code>phiTaylorPhy()</code>)
%> 
%> @param  p   The polynomial degree
%> @retval ret A matrix of two-dimensional multi-indices.
%>             @f$[\frac{p(p+1)}{2} \times 2]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%> 
%> @par License
%> @parblock
%> This program is free software: you can redistribute it and/or modify
%> it under the terms of the GNU General Public License as published by
%> the Free Software Foundation, either version 3 of the License, or
%> (at your option) any later version.
%>
%> This program is distributed in the hope that it will be useful,
%> but WITHOUT ANY WARRANTY; without even the implied warranty of
%> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%> GNU General Public License for more details.
%>
%> You should have received a copy of the GNU General Public License
%> along with this program.  If not, see <http://www.gnu.org/licenses/>.
%> @endparblock
%
function mult = multiindex(p)
mult = zeros(p * (p+1) / 2, 2);
mult(1,:) = [0, 0];
for ord = 1 : p
  offset = ord * (ord+1) / 2;
  for i = 1 : ord + 1
    mult(offset + i, :) = mult(1, :) + [ord - i + 1, i - 1];
  end
end
end