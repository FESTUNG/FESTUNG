% Evaluates the i-th one-dimensional basis function.

%===============================================================================
%> @file phi1D.m
%>
%> @brief Evaluates the i-th basis function.
%===============================================================================
%>
%> @brief Evaluates the @f$i@f$-th basis function on the reference interval 
%>        @f$[0,1]@f$ at points specified by a list of @f$\hat{x}@f$.
%>
%> @param  i   The index of the basis function.
%> @param  X   A list of @f$\hat{x}@f$ coordinates.
%> @retval ret The @f$i@f$-th basis function in all points specified by
%>             <code>X</code>. It holds <code>size(X) == size(ret)</code>.
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
function ret = phi1D(i, X)
switch i
  case 1,  ret = ones(size(X));
  case 2,  ret = sqrt(3) * (1 - 2 * X);
  case 3,  ret = sqrt(5) * ( (6 * X - 6) .* X + 1 );
  case 4,  ret = sqrt(7) * ( ( ( 20 * X - 30 ) .* X + 12 ) .* X - 1 );
end % switch
end % function