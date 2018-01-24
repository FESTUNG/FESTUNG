% Assembles matrices containing a two-dimensional basis function evaluated
% in all quadrature points of the edges.

%===============================================================================
%> @file assembleMatEdgePhiIntPerQuad.m
%>
%> @brief Assembles matrices containing a two-dimensional basis function 
%>        evaluated in all quadrature points of the edges.
%===============================================================================
%>
%> @brief Assembles matrices containing a two-dimensional basis function 
%>        evaluated in all quadrature points of the edges.
%>
%> The matrices @f$\mathsf{{S}}^n \in \mathbb{R}^{KN\times KR}@f$are defined as 
%> @f[
%> [\mathsf{{S}}^n]_{(k-1)N+i,(k-1)N+r} =
%>  w_r \varphi_{ki}(\mathbf{q}_r) \,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>domainRectTrap()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times n_\mathrm{edges}]@f$
%> @param refEdgePhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{{S}}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgeTetraPhiIntPerQuad()</code>.
%>                    @f$[N \times R \times n_\mathrm{edges}]@f$
%> @retval ret        The assembled matrices @f$[ n_\mathrm{edges}\times1 \mathrm{cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2017
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
function ret = assembleMatEdgePhiIntPerQuad(g, markE0T, refEdgePhiIntPerQuad)
K = g.numT;  nEdges = size(g.E0T, 2);
ret = cell(nEdges, 1);
for n = 1 : nEdges
	ret{n} = kron( spdiags(markE0T(:, n) .* g.areaE0T(:, n), 0, K, K), refEdgePhiIntPerQuad(:, :, n) );
end % for nn
end % function
