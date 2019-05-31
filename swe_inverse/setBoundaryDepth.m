% Set the boundary values for the bathymetry.

%===============================================================================
%> @file
%>
%> @brief Set the boundary values for the bathymetry.
%===============================================================================
%>
%> @brief Set the boundary values for the bathymetry.
%>
%> @param  pd           A struct with problem parameters and precomputed
%>                      fields, as provided by configureProblem() and 
%>                      preprocessProblem(). @f$[\text{struct}]@f$
%>
%> @retval ret					The discrete data representing the bathymetry on
%>											Dirichlet type boundaries. 
%>                      @f$[3 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function ret = setBoundaryDepth(pd)
K = pd.K;
[Q, ~] = quadRule1D(2*pd.p+1); numQuad1D = length(Q);
ret = cell(3,1);
for nn = 1 : 3
  ret{nn} = reshape(pd.basesOnQuad.phi1D{2*pd.p+1}(:,:,nn) * pd.zbExact(:,:).', numQuad1D * K, 1);
end % for
end % function
