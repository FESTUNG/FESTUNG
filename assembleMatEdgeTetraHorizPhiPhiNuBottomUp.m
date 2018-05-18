% Assembles a matrix containing integrals over interior edges of products of
% two basis functions with the x2-component of the edge normal, where one of the 
% basis functions always comes from the element below the edge.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals over interior edges of 
%>        products of two basis functions with the @f$x^2@f$-component of the
%>        edge normal, where one of the basis functions always stems from the 
%>        element below the edge.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals over interior edges of 
%>        products of two basis functions with the @f$x^2@f$-component of the
%>        edge normal, where one of the basis functions always stems from the 
%>        element below the edge.
%>
%> The entries of the matrix have essentially the same form as for
%> core/assembleMatEdgePhiPhiNu.m with the difference that
%> - only horizontal edges are considered (i.e., @f$n=1,2@f$)
%> - @f$\varphi_j@f$ always stems from the element below the current edge
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridDataTetra()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    <code>logical</code> arrays that mark each elements
%>                    edges on which the matrix blocks should be
%>                    assembled @f$[K \times 4]@f$
%> @param refEdgePhiIntPhiInt Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiInt()</code>.
%>                    @f$[N \times N \times 4]@f$
%> @param refEdgePhiIntPhiExt Local matrix 
%>                    @f$\hat{\mathsf{Q}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiExt()</code>.
%>                    @f$[N \times N \times 4]@f$
%>
%> @retval ret        The assembled matrix @f$[KN\times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018.
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
function ret = assembleMatEdgeTetraHorizPhiPhiNuBottomUp(g, markE0T, refEdgePhiIntPhiInt, refEdgePhiIntPhiExt)
areaNuMarkE0T = { markE0T(:,1) .* g.areaE0T(:,1) .* g.nuE0T(:,1,2), ...
                  markE0T(:,2) .* g.areaE0T(:,2) .* g.nuE0T(:,2,2) };
ret = kron(spdiags(areaNuMarkE0T{2}, 0, g.numT, g.numT ), refEdgePhiIntPhiInt(:,:,2)) + ...
      kron(bsxfun(@times, g.markE0TE0T{1}, areaNuMarkE0T{1}), refEdgePhiIntPhiExt(:,:,1));
end  % function