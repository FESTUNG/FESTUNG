% Obtain local edge index np in neighboring element from local edge index 
% nn for quadrilateral elements.

%===============================================================================
%> @file
%>
%> @brief Obtain local edge index np in neighboring element from local edge
%>        index nn for quadrilateral elements.
%===============================================================================
%>
%> @brief Obtain local edge index np in neighboring element from local edge
%>        index nn for quadrilateral elements.
%>
%> Due to the structured nature of quadrilateral meshes, the local edge
%> index @f$n^+@f$ in the neighbouring element is implicitly given by the
%> local edge index @f$n^-@f$.
%>
%> Since the direction of integration does not change in our quadrilateral
%> meshes there is no need for a mapping @f$\theta@f$ (theta()) as for
%> triangular meshes. However, there may be a need to identify the local
%> edge index for an edge in the neighboring element.
%>
%> This function can be used to compute one from the other.
%> 
%> @param  nn    Local edge index @f$n^-@f$.
%>
%> @retval np    Local edge index @f$n^+@f$.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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
function np = mapLocalEdgeTetra(nn)
mapE0E = [2 1 4 3];
np = mapE0E(nn);
end % function