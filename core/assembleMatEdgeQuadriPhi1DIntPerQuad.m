% Assembles matrices containing a one-dimensional basis function evaluated
% in all quadrature points of the edges.

%===============================================================================
%> @file
%>
%> @brief Assembles matrices containing a one-dimensional basis function 
%>        evaluated in all quadrature points of the edges.
%===============================================================================
%>
%> @brief Assembles matrices containing a one-dimensional basis function 
%>        evaluated in all quadrature points of the edges.
%>
%> The matrices @f$\mathsf{{S}}^n \in \mathbb{R}^{KN\times \overline{K}R}@f$
%> are defined as 
%> @f[
%> [\mathsf{{S}}^n]_{(k-1)N+i,(k-1)N+r} =
%>  w_r \phi_{ki}(q_r) \,.
%> @f]
%>
%> @param  g2D        The lists describing the geometric and topological 
%>                    properties of a quadrilateral triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  g1D        The lists describing the geometric and topological 
%>                    properties of a matching one-dimensional triangulation
%>                    (see <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhi1DIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeQuadriPhi1DIntPerQuad()</code>.
%>                    @f$[N \times R \times 4]@f$
%> @retval ret        The assembled matrices @f$[4\times1 \mathrm{cell}]@f$
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
function ret = assembleMatEdgeQuadriPhi1DIntPerQuad(g2D, g1D, markE0T, refEdgePhi1DIntPerQuad)
K = g2D.numT; barK = g1D.numT;
ret = cell(4,1);
for nn = 1 : 4
	ret{nn} = kron( spdiags(g1D.markT2DT.' * ( markE0T(:, nn) .* g2D.areaE0T(:, nn) ), 0, barK, K), ...
                  refEdgePhi1DIntPerQuad(:, :, nn) );
end % for nn
end % function
