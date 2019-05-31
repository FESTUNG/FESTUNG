% Computes the two indices of the underlying one-dimensional basis
% functions for a given linearized index of a two-dimensional basis tensor
% product basis function.

%===============================================================================
%> @file
%>
%> @brief Computes the two indices of the underlying one-dimensional basis
%>        functions for a given linearized index of a two-dimensional basis tensor
%>        product basis function.
%===============================================================================
%>
%> @brief Computes the two indices of the underlying one-dimensional basis
%>        functions for a given linearized index of a two-dimensional basis tensor
%>        product basis function.
%>
%> This is the inverse of <code>mapTensorProductIndex()</code>.
%> 
%> @param  i     Linearized index of the two-dimensional basis function.
%> @retval m,n   Index of one-dimensional basis functions.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function [m, n] = mapTensorProductIndexInv(i)
p = ceil(sqrt(i));
local_ind = p-(i-(p-1).^2);
m = p + min(local_ind, 0);
n = p - max(local_ind, 0);
end

