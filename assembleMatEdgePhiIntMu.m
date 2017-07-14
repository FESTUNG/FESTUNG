% TODO
% Assembles a matrix containing integrals of products of an interior basis
% with an edge basis function.

%===============================================================================
%> @file assembleMatEdgePhiIntMu.m
%>
%> @brief Assembles a matrix containing integrals of products of an interior
%>        basis with an edge basis function.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals of products of an interior
%>        basis with an edge basis function.
%>
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param refEdgePhiIntMu  Local matrix @f$\hat{\mathsf{R}}_{\mu}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntMu()</code>.
%>                    @f$[N \times \bar{N}]@f$
%> @retval ret        The assembled matrix @f[KN \times \bar{K}\bar{N}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
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
function ret = assembleMatEdgePhiIntMu(g, markE0T, refEdgePhiIntMu)
K = g.numT;  Kedge = g.numE; 
[N, Nmu, ~, ~] = size(refEdgePhiIntMu);

% Check function arguments that are directly used
validateattributes(markE0T, {'logical'}, {'size', [K 3]}, mfilename, 'markE0T');
validateattributes(refEdgePhiIntMu, {'numeric'}, {'size', [N, Nmu, 3, 2]}, mfilename, 'refEdgePhiIntMu');

ret = sparse(K * N, Kedge * Nmu);
for n = 1 : 3
  for l = 1 : 2
    Rkn = markE0T(:,n) .*  g.markSideE0T(:, n, l) .* g.areaE0T(:, n);
    ret = ret + kron(sparse(1 : K, g.E0T(:, n), Rkn, K, Kedge), refEdgePhiIntMu(:, :, n, l));
  end
end % for
end % function
