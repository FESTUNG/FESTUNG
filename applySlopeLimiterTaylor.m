% Applies a slope limiter to a discrete function given in Taylor basis.
%
%===============================================================================
%> @file applySlopeLimiterTaylor.m
%>
%> @brief Applies a slope limiter to a discrete function given in Taylor basis.
%===============================================================================
%>
%> @brief Applies a slope limiter to a discrete function given in Taylor basis.
%>
%> It applies the slope limiting operator @f$\mathsf{\Phi}^\mathrm{Taylor}@f$
%> to a discrete function @f$c_h: \Omega \rightarrow \mathbb{R}@f$ given as
%> representation matrix, e.g., obtained by projection into Taylor basis
%> (cf. <code>projectDataDisc2DataTaylor()</code>).
%>
%> The type of slope limiter is determined by the parameter
%> <code>type</code>, which can be one of the following:
%> 
%> - <code>'linear'</code>: See <code>applySlopeLimiterTaylorLinear()</code>
%> - <code>'hierarch_vert'</code>: See <code>applySlopeLimiterTaylorHierarchicalVertex()</code>
%> - <code>'strict'</code>: See <code>applySlopeLimiterTaylorStrict()</code>
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataTaylor The representation matrix of the unlimited function 
%>                    @f$c_h@f$. @f$[K \times N]@f$
%> @param  markV0TbdrD <code>logical</code> arrays that mark each triangles
%>                    (Dirichlet boundary) vertices on which additional
%>                    function values should be considered during the 
%>                    slope limiting routine. @f$[K \times 3]@f$
%> @param  dataV0T    The function values for (Dirichlet boundary) vertices
%>                    specified by <code>markV0TbdrD</code>. @f$[K \times 3]@f$
%> @param  type       The type of slope limiter to be used. [<code>string</code>]
%> @retval dataTaylorLim   The representation matrix of the limited function
%>                    @f$\mathsf{\Phi}^\mathrm{Taylor}c_h@f$. @f$[K \times N]@f$
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
function dataTaylorLim = applySlopeLimiterTaylor(g, dataTaylor, markV0TbdrD, dataV0T, type)
% Set default limiter
if nargin < 5
  type = 'linear';
end
% Apply chosen limiter
switch type
  case 'linear'
    dataTaylorLim = applySlopeLimiterTaylorLinear(g, dataTaylor, markV0TbdrD, dataV0T);
  case 'hierarch_vert'
    dataTaylorLim = applySlopeLimiterTaylorHierarchicalVertex(g, dataTaylor, markV0TbdrD, dataV0T);
  case 'strict'
    dataTaylorLim = applySlopeLimiterTaylorStrict(g, dataTaylor, markV0TbdrD, dataV0T);
  otherwise
    error('Unknown limiter type');
end
end

