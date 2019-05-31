% Evaluates the i-th linear Lagrangian basis function.
%
%===============================================================================
%> @file phi.m
%>
%> @brief Evaluates the i-th linear Lagrangian basis function.
%===============================================================================
%>
%> @brief Evaluates the @f$i@f$-th linear Lagrangian basis function on the 
%>        reference triangle @f$\hat{T}@f$ at points specified by a list of
%>        @f$\hat{x}^1@f$ and @f$\hat{x}^2@f$ coordinates .
%>
%> @param  i   The index of the basis function.
%> @param  X1  A list of @f$\hat{x}^1@f$ coordinates.
%> @param  X2  A list of @f$\hat{x}^2@f$ coordinates.
%> @retval ret The @f$i@f$-th basis function in all points specified by
%>             <code>X1</code>, <code>X2</code>. It holds <code>size(X1) ==
%>             size(X2) == size(ret)</code>.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = phiLagr(i, X1, X2)
switch i
  case 1, ret = 1 - X1 - X2;
  case 2, ret = X1;
  case 3, ret = X2;
end % switch
end % function
