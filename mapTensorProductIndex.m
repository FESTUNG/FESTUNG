% Computes the one-dimensional (linearized) index of a tensor product basis 
% for two given indices of the underlying basis functions.

%===============================================================================
%> @file mapTensorProductIndex.m
%>
%> @brief Computes the one-dimensional (linearized) index of a tensor product
%>         basis for two given indices of the underlying basis functions.
%===============================================================================
%>
%> @brief Computes the one-dimensional (linearized) index of a two-dimensional
%>        tensor product basis function for two given indices of the associated 
%>        one-dimensional basis functions.
%>
%> This gives the linearized index of a basis function, as provided by
%> <code>phiTensorProduct()</code>, for given indices of the underlying 
%> basis functions.
%>
%> The index is chosen such that basis functions of a polynomial degree
%> appear consecutively.
%>
%> The mapping up to polynomial degree three is the following:
%> <table>
%> <tr><th> p </th> <td> 0 </td> <td colspan="3"> 1 </td> <td colspan="5"> 2 </td> 
%>   <td colspan="7"> 3 </td></tr>
%> <tr><th> [m,n] </th> <td> [1,1] </td> <td> [2,1] </td> <td> [2,2] </td> <td> [1,2] </td> 
%>   <td> [3,1] </td> <td> [3,2] </td> <td> [3,3] </td> <td> [2,3] </td> <td> [1,3] </td>
%>   <td> [4,1] </td> <td> [4,2] </td> <td> [4,3] </td> <td> [4,4] </td> <td> [3,4] </td>
%>   <td> [2,4] </td> <td> [1,4] </td></tr>
%> <tr><th> i </th> <td> 1 </td> <td> 2 </td> <td> 3 </td> <td> 4 </td> 
%>   <td> 5 </td> <td> 6 </td> <td> 7 </td> <td> 8 </td> <td> 9 </td> 
%>   <td> 10 </td> <td> 11 </td> <td> 12 </td> <td> 13 </td> <td> 14 </td> 
%>   <td> 15 </td> <td> 16 </td> </tr>
%> </table>
%> 
%> @param  m, n  Index of one-dimensional basis functions.
%> @retval ret   Linearized index of the two-dimensional basis function.
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
function ret = mapTensorProductIndex(m, n)
p = max(m,n);
ret = (p - 1) .* (p - 1) + p - m + n;
end

