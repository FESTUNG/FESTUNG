% Assembles two matrices, each containing integrals of products of the continuous advection velocity, a basis function and a (spatial) derivative of a basis function.

%===============================================================================
%> @file
%>
%> @brief Assembles two matrices, each containing integrals of products of the advection velocity, a basis function and a (spatial) derivative of a basis function.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{G}^m, m \in \{1,2\}@f$ containing integrals of products of the advection velocity, a basis function and a (spatial) derivative of a basis function.
%>
%> @param  g             		The lists describing the geometric and topological 
%>		                      properties of a triangulation (see 
%>      		                <code>generateGridData()</code>) 
%>      		                @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhiPerQuad Precomputed contributions to the local matrix 
%> 						            	@f$\hat{\mathsf{G}}@f$ 
%>                          as provided by <code>integrateRefElemDphiPhiPerQuad()</code>.
%>                          @f$[2 \times 1 \text{ cell}]@f$
%> @param funcCont1        A function handle for the continuous velocity in x-direction.
%> @param funcCont2        A function handle for the continuous velocity in y-direction.
%> @param qOrd		         The order of the quadrature rule.
%>                    
%>                    
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%> 
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
function ret = assembleMatElemDphiDphiFuncCont(g, refElemDphiDphiPerQuad, funcCont, qOrd)
K = g.numT;
[N, ~, R] = size(refElemDphiDphiPerQuad{1,1});
validateattributes(refElemDphiDphiPerQuad, {'cell'}, {'size', [2 2]}, mfilename, 'refElemDphiDphiPerQuad');
validateattributes(refElemDphiDphiPerQuad{1,1}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiDphiPerQuad{1,1}');
validateattributes(refElemDphiDphiPerQuad{1,2}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiDphiPerQuad{1,2}');
validateattributes(refElemDphiDphiPerQuad{2,2}, {'numeric'}, {'size', [N N R]}, mfilename, 'refElemDphiDphiPerQuad{2,2}');
validateattributes(funcCont, {'function_handle'}, {'size', [1 1]}, mfilename, 'funcCont');

if nargin < 4, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1); end
[Q1, Q2, ~] = quadRule2D(qOrd);

% Assemble matrix
ret = { zeros(K*N, N), zeros(K*N, N) };
for r = 1 : R
  valOnQuad = funcCont(g.mapRef2Phy(1, Q1(r), Q2(r)), g.mapRef2Phy(2, Q1(r), Q2(r)));
  ret{1} = ret{1} + kron(g.B(:,2,2) .* g.B(:,2,2) .* valOnQuad, refElemDphiDphiPerQuad{1,1}(:, :, r)) ...
                  - kron(2 * g.B(:,2,2) .* g.B(:,2,1) .* valOnQuad, refElemDphiDphiPerQuad{1,2}(:, :, r)) ...
                  + kron(g.B(:,2,1) .* g.B(:,2,1) .* valOnQuad, refElemDphiDphiPerQuad{2,2}(:, :, r));

  ret{2} = ret{2} + kron(g.B(:,1,2) .* g.B(:,1,2) .* valOnQuad, refElemDphiDphiPerQuad{1,1}(:, :, r)) ...
                  - kron(2 * g.B(:,1,2) .* g.B(:,1,1) .* valOnQuad, refElemDphiDphiPerQuad{1,2}(:, :, r)) ...
                  + kron(g.B(:,1,1) .* g.B(:,1,1) .* valOnQuad, refElemDphiDphiPerQuad{2,2}(:, :, r));
end % for
ret{1} = kronVec(spdiags(1 ./ g.detJ0T, 0, K, K), ret{1});
ret{2} = kronVec(spdiags(1 ./ g.detJ0T, 0, K, K), ret{2});
end % function
