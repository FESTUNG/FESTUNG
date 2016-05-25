% Compute the linear index corresponding to a two-dimensional multiindex.
%
%===============================================================================
%> @file mult2ind.m
%>
%> @brief Compute the linear index corresponding to a two-dimensional multiindex.
%===============================================================================
%>
%> @brief Compute the linear index corresponding to a two-dimensional multiindex.
%>
%> For a given two-dimensional multiindex @f$\mathbf{a} = [a^1, a^2]^T@f$
%> the corresponding linear index is computed, using the following ordering:
%>
%>        a^T   |  ind
%>     ---------+-------
%>       [0,0]  |   1
%>       [1,0]  |   2
%>       [0,1]  |   3
%>       [2,0]  |   4
%>       [1,1]  |   5
%>       [0,2]  |   6
%>       [3,0]  |   7
%>       [2,1]  |   8
%>       [1,2]  |   9
%>       [0,3]  |  10
%>             ...
%>
%> This is useful, e.g., to find the index of a certain derivative in the
%> Taylor basis functions.
%>
%> @param a      Two-dimensional multi-index as a row vector. Multiple
%>               multi-indices can be given as multiple rows. 
%>               @f$[n \times 2]@f$
%> 
%> @return The linear index. @f$[n \times 1]@f$
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
function ind = mult2ind(a)
assert(size(a, 2) == 2);
p = sum(a, 2);
N = p .* (p + 1) / 2;
ind = N + 1 + a(:, 2);
end