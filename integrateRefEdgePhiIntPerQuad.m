% Evaluates a basis function in all quadrature points of all edges of the 
% reference square.

%===============================================================================
%> @file integrateRefEdgeTetraPhiIntPerQuad.m
%>
%> @brief Evaluates a basis function in all quadrature points of all edges of 
%>        the reference square.
%===============================================================================
%>
%> @brief Evaluates a basis function in all quadrature points of all edges of 
%>        the reference square.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{R}}^\mathrm{diag} 
%>    \in \mathbb{R}^{N\times R\times n_\mathrm{edges}}@f$, which is defined by
%> @f[
%> \left[\hat{\mathsf{R}}^\mathrm{diag}\right]_{i,r,n} \;:=\;
%> \hat{\varphi}_i\circ\hat{\mathbf{\gamma}}_n(q_r)
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in <code>gammaMap()</code> or
%> <code>gammaMapTetra()</code> and the quadrature points @f$q_r@f$ given by
%> <code>quadRule1D()</code>
%>
%> @param  N            The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @param  qOrd         The order of the quadrature rule to be used
%> @retval ret  The computed array @f$[N\times R\times n_\mathrm{edges}]@f$
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
function ret = integrateRefEdgePhiIntPerQuad(N, basesOnQuad, qOrd)
[~, W] = quadRule1D(qOrd); R = length(W);
nEdges = size(basesOnQuad.phi1D{qOrd}, 3);
ret = zeros(N, R, nEdges);
for n = 1 : nEdges
  for i = 1 : N
    ret(i, :, n) = W .* basesOnQuad.phi1D{qOrd}(:, i, n).';
  end  % for i
end  % for n
end  % function