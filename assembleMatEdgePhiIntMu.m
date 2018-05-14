% Assembles a matrix containing integrals of products of an interior
% basis with an edge basis function.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals of products of an interior
%> basis with an edge basis function.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals of products of an interior
%> basis with an edge basis function.
%>
%> The matrix @f$\mathsf{R}_{\mu} \in \mathbb{R}^{KN \times \bar{K} \bar{N}}@f$
%> with @f$K@f$ the number of elements, @f$N@f$ the number of local degrees of 
%> freedom on the element, @f$\bar{K}@f$ the number of edges and @f$\bar{N}@f$ 
%> the number of local degrees of freedom on an edge is given componentwise by
%> 
%> @f[
%>   [\mathsf{R}_{\mu}]_{(k-1)N+i,(\bar{k}-1)\bar{N}+j} = \sum_{E_{kn} \in \partial{{T_{k}}} \cap \ensuremath{\mathcal{E}}_{\text{int}} } \int_{E_{kn}}  \mu_{knj} \, \varphi_{ki} \, \text{d}s.
%> @f]
%> 
%> The assembled matrix is given as
%> 
%> @f[
%> \mathsf{R}_{\mu} 
%> = \sum_{n=1}^3 \sum_{l=1}^2 \left( \begin{bmatrix}
%>   \ensuremath{|E_{1n}|}\, \delta_{E_{1n}\in\ensuremath{\mathcal{E}}_{\text{int}}} & & \\
%>   & \ddots & \\
%>   & & \ensuremath{|E_{Kn}|}\, \delta_{E_{Kn}\in\ensuremath{\mathcal{E}}_{\text{int}}}
%> \end{bmatrix} \, \Delta_n \right) \otimes [\mathsf{\hat{R}}_{\mu}]_{:,:,n,l} 
%> @f]
%> 
%> where @f$\delta_{E_{kn}\in\ensuremath{\mathcal{E}}_{\text{int}}}@f$ (\p markE0T) denotes
%> a Kronecker delta, @f$\Delta_n@f$ is the permutation matrix described in 
%> <code>assembleMatEdgePhiIntMuVal()</code> and @f$\otimes@f$ denotes the
%> Kronecker product. The local matrices 
%> @f[
%> [\mathsf{\hat{R}}_{\mu}]_{i,j,n,l} := \int_{0}^{1} \hat{\varphi}_{i} \circ \boldsymbol{\hat{\gamma}}_{n}(s) \, \hat{\mu}_{j} \circ \hat{\beta}_{kn}(s) \, \text{d}s 
%> @f]
%> are precomputed for all combinations of element and edge test functions and
%> given to the routine through \p refEdgePhiIntMu.
%> 
%> The matrix @f$\mathsf{T}@f$ can be constructed by this routine as
%> @f[
%> \mathsf{R}_{\mu} = {\mathsf{T}}^\mathrm{T}
%> @f]
%> holds.
%> 
% %> @f[
% %> \int_{E_{kn}}  \varphi_{ki} \, \mu_{knj} \, \text{d}s	=
% %> \ensuremath{|E_{kn}|}
% %> \underbrace{\int_{0}^{1} \hat{\varphi}_{i} \circ \boldsymbol{\hat{\gamma}}_{n}(s) \, \hat{\mu}_{j} \circ \hat{\beta}_{kn}(s) \, \text{d}s}_{\eqqcolon [\mathsf{\hat{R}}_{\mu}]_{i,j,n,l}}\,,
% %> @f]
% %> 
% %> @f[
% %>   \mathsf{R}_{\mu,E_{kn}} = \int_{E_{kn}} 
% %>   \begin{bmatrix}
% %>     \mu_{kn1}\,\varphi_{k1} & \ldots & \mu_{kn\bar{N}}\,\varphi_{k1} \\
% % %>     \vdots           & \ddots & \vdots \\
% %>     \mu_{kn1}\,\varphi_{kN} & \ldots & \mu_{kn\bar{N}}\,\varphi_{kN} 
% %>   \end{bmatrix}
% %>   \text{d}s.
% %> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param refEdgePhiIntMu  Local matrices @f$\hat{\mathsf{R}}_{\mu}@f$ as provided
%>                    by <code>integrateRefEdgePhiIntMu()</code>.
%>                    @f$[N \times \bar{N} \times 3 \times 2]@f$
%> @retval ret        The assembled matrix @f$[KN \times \bar{K}\bar{N}]@f$
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
  end % for
end % for
end % function
