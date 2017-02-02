% Evaluates the i-th basis function.
%
%===============================================================================
%> @file phi.m
%>
%> @brief Evaluates the i-th basis function.
%===============================================================================
%>
%> @brief Evaluates the @f$i@f$-th basis function on the reference triangle
%>        @f$\hat{T}@f$ at points specified by a list of @f$\hat{x}^1@f$
%>        and @f$\hat{x}^2@f$ coordinates .
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
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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

%Mapping basis function from [-1,1] onto [0,1]
xi = 2 * ( X ) - 1;

switch i
    case 1,  ret = 1.;
    case 2,  ret = xi;
    case 3,  ret = 0.5 * (3.*xi.^2 - 1);
    case 4,  ret = 0.125 * ( 35 * xi.^4 - 30*xi.^2 + 3);
    case 5,  ret = 0.125 * ( 63 * xi.^5 - 70*xi.^3 + 15*xi);
    case 6,  ret = 1/16 * ( 231 * xi.^6 - 315*xi.^4 + 105*xi.^2 - 5);
end % switch
end % function
