% Assembles a matrix containing integrals of products of two edge basis
% functions.

%===============================================================================
%> @file assembleMatEdgeMuMu.m
%>
%> @brief Assembles a matrix containing integrals of products of two edge basis 
%>        functions. This corresponds to a mass matrix.
%===============================================================================
%>
%> @brief Assembles a mass matrix @f$\mathsf{M}_{\mu}@f$
%>        containing integrals of products of two edge basis functions.
%>
%> The matrix @f$\mathsf{M}_\mu \in \mathbb{R}^{\bar{K}\bar{N}\times \bar{K}\bar{N}}@f$,
%> with @f$\bar{K}@f$ the number of edges and @f$\bar{N}@f$ the number of local degrees
%> of freedom on an edge, is block diagonal and defined component-wise by
%> @f[
%>   [\mathsf{M}_mu]_{(\bar{k}-1)\bar{N}+i,(\bar{k}-1)\bar{N}+j} = 
%>      \sum_{E_{kn}\in\mathcal{E}} \int_{E_{kn}} \mu_{kni} \mu_{knj} \mathrm{d}s \,,
%> @f]
%> where @f$\mathcal{E}@f$ is the set of edges to be considered.
%>
%> The integral is backtransformed to the reference interval @f$[0,1]@f$, defined as
%> @f[
%>   \int_{E_{kn}} \varphi_{ki} \mu_{knj} \mathrm{ds} = ...
%> @f]
%>
%> TODO (see assembleMatEdgePhiPhi).
%> 
%> All other entries are zero.
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param refElemMuMu Local matrix @f$\hat{\mathsf{M}}_{\mu}@f$ as provided
%>                    by <code>integrateRefElemMuMu()</code>.
%>                    @f$[\bar{N} \times \bar{N}]@f$
%> @retval ret        The assembled matrix @f$[\bar{K}\bar{N} \times \bar{K}\bar{N}]@f$
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
function ret = assembleMatEdgeMuMu(g, markE0T, refEdgeMuMu)
Nmu = size(refEdgeMuMu, 1);
Kedge = g.numE;
validateattributes(refEdgeMuMu, {'numeric'}, {'size', [Nmu Nmu]}, mfilename, 'refEdgeMuMu');
ret = sparse(Kedge * Nmu, Kedge * Nmu);
for n = 1 : 3
  Kkn = g.areaE0T(:, n) .*  markE0T(:, n) ;
  ret = ret + kron(sparse(g.E0T(:, n), g.E0T(:, n), Kkn, Kedge, Kedge), refEdgeMuMu);
end
end % function
