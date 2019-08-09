% Assembles a matrix, containing integrals of products of two basis
% functions and a continuous function.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix, containing integrals of products of two basis 
%>        functions and a continuous function.
%===============================================================================
%>
%> @brief Assembles matrix @f$\mathsf{A}@f$
%>        containing integrals of products of two basis functions and a
%>        continuous function.
%>
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemPhiPhi Local matrix @f$\hat{\mathsf{M}}@f$ as provided
%>                    by <code>integrateRefElemPhiPhi()</code>.
%>                    @f$[N \times N]@f$
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
function ret = assembleMatElemPhiPhiFuncCont(g, refElemPhiPhiPerQuad, funcCont, qOrd)
K = g.numT;
[N, ~, R] = size(refElemPhiPhiPerQuad);

validateattributes(refElemPhiPhiPerQuad, {'numeric'}, {'size', [N, N, R]}, mfilename, 'refElemPhiPhiPerQuad');
validateattributes(funcCont, {'function_handle'}, {'size', [1 1]}, mfilename, 'funcCont');

if nargin < 4, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); end
[Q1, Q2, ~] = quadRule2D(qOrd);

% Assemble matrix
ret = sparse(K*N, K*N);
for r = 1 : R
  valOnQuad = funcCont(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
  ret = ret + kron(spdiags(g.detJ0T .* valOnQuad, 0, K, K), refElemPhiPhiPerQuad(:, :, r));
end % for
end % function
