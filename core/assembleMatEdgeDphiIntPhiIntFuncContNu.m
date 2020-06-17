% Assembles a matrix containing integrals over interior edges of products of
% two basis functions with the upwind value of a, for each quadrature point 
% specified, function.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals over interior edges of products of
%>        two basis functions with the upwind value of a, for each quadrature 
%>        point specified, function.
%===============================================================================
%>
%> @brief Assembles matrix @f$\mathsf{R}@f$ containing integrals over
%>        edges of products of two basis functions with the upwind value of 
%>        a, for each quadrature point specified, function.
%>
%> @param g           The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param markE0Tbdr  <code>logical</code> arrays that mark each triangles
%>                    (boundary) edges on which the matrix blocks should be
%>                    assembled @f$[K \times 3]@f$
%> @param refEdgePhiIntPhiIntPerQuad  Local matrix 
%>                    @f$\hat{\mathsf{R}}^\text{diag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiIntPerQuad()</code>.
%>                    @f$[N \times N \times 3 \times R]@f$
%> @param refEdgePhiIntPhiExtPerQuad Local matrix 
%>                    @f$\hat{\mathsf{R}}^\text{offdiag}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntPhiExtPerQuad()</code>.
%>                    @f$[N \times N \times 3 \times 3 \times R]@f$
%> @param funcCont    A function handle for the continuous function
%> @param qOrd        (optional) Order of quadrature rule to be used.
%                     Defaults to @f$2p+1@f$.
%>
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2019 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>                      
%> @author Balthasar Reuter, 2019
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
function ret = assembleMatEdgeDphiIntPhiIntFuncContNu(g, markE0Tbdr, refEdgeDphiIntPhiIntPerQuad, funcCont, qOrd)
% Extract dimensions
K = g.numT;
[N, ~, ~, R] = size(refEdgeDphiIntPhiIntPerQuad{1});
if nargin < 5, p = (sqrt(8*N+1)-3)/2; qOrd = 2*p + 1; end
[Q, ~] = quadRule1D(qOrd);

% Check function arguments that are directly used
validateattributes(markE0Tbdr, {'logical'}, {'size', [K 3]}, mfilename, 'markE0Tbdr');
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(refEdgeDphiIntPhiIntPerQuad, {'cell'}, {'numel', 2});
validateattributes(refEdgeDphiIntPhiIntPerQuad{1}, {'numeric'}, {'size', [N N 3 R]});
validateattributes(refEdgeDphiIntPhiIntPerQuad{2}, {'numeric'}, {'size', [N N 3 R]});

% Assemble matrices
ret = { sparse(K*N, K*N), sparse(K*N, K*N) };
for nn = 1 : 3
  [Q1, Q2] = gammaMap(nn, Q);
  valOnQuad = funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2));
  
  Rkn = markE0Tbdr(:, nn) .* g.areaE0T(:, nn) ./ g.detJ0T;
  % Diagonal blocks
  for r = 1 : R
    markRknValOnQuad = Rkn .* valOnQuad(:, r);
    ret{1} = ret{1} + kron(spdiags(markRknValOnQuad .* g.nuE0T(:, nn, 1) .* g.B(:, 2, 2), 0, K, K), ...
                           refEdgeDphiIntPhiIntPerQuad{1}(:, :, nn, r)) ...
                    - kron(spdiags(markRknValOnQuad .* g.nuE0T(:, nn, 1) .* g.B(:, 2, 1), 0, K, K), ...
                           refEdgeDphiIntPhiIntPerQuad{2}(:, :, nn, r));
    ret{2} = ret{2} - kron(spdiags(markRknValOnQuad .* g.nuE0T(:, nn, 2) .* g.B(:, 1, 2), 0, K, K), ...
                           refEdgeDphiIntPhiIntPerQuad{1}(:, :, nn, r)) ...
                    + kron(spdiags(markRknValOnQuad .* g.nuE0T(:, nn, 2) .* g.B(:, 1, 1), 0, K, K), ...
                           refEdgeDphiIntPhiIntPerQuad{2}(:, :, nn, r));
  end
end % for
end % function
