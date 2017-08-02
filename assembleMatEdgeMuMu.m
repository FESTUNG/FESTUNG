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
%>   [\mathsf{M}_{\mu}]_{(\bar{k}-1)\bar{N}+i,(\bar{k}-1)\bar{N}+j} = 
%>      \sum_{E_{kn}\in\mathcal{E}} \int_{E_{kn}} \mu_{kni} \mu_{knj} \mathrm{d}s \,,
%> @f]
%> where @f$\mathcal{E}@f$ is the set of edges to be considered.
%>
%> Two different mass matrices have to be considered. A mass matrix resulting 
%> from edges from the interior of the mesh
%> @f$E_{kn} \in \ensuremath{\mathcal{E}}_{\text{int}}@f$ 
%> @f[
%> \mathsf{\bar{M}}_{\mu} = 
%> \sum_{n=1}^3 \left( {(\mathsf{\Delta}_n)}^\mathrm{T} \, \begin{bmatrix}
%> \ensuremath{|E_{1n}|}\, \delta_{E_{1n}\in\ensuremath{\mathcal{E}}_{\text{int}}} & & \\
%>  & \ddots &  \\
%>   &  & \ensuremath{|E_{Kn}|} \,\delta_{E_{Kn}\in\ensuremath{\mathcal{E}}_{\text{int}}}
%> \end{bmatrix} \, \mathsf{\Delta}_n \right) \otimes  \hat{\mathsf{M}}_\mu
%> @f]
%> and edges representing the domain boundary @f$E_{kn} \in \ensuremath{\mathcal{E}}_{\text{bc}}@f$ 
%> @f[
%> \tilde{\mathsf{M}}_{\mu} =
%> \sum_{n=1}^3 \left( {(\mathsf{\Delta}_n)}^\mathrm{T} \, \begin{bmatrix}
%> \ensuremath{|E_{1n}|}\, \delta_{E_{1n}\in\ensuremath{\mathcal{E}}_{\text{bc}}} & & \\
%>  & \ddots &  \\
%>   &  & \ensuremath{|E_{Kn}|} \,\delta_{E_{Kn}\in\ensuremath{\mathcal{E}}_{\text{bc}}}
%> \end{bmatrix} \, \mathsf{\Delta}_n \right) \otimes  \hat{\mathsf{M}}_\mu
%> @f]
%> where @f$\delta_{E_{1n}\in\ensuremath{\mathcal{E}}_{\text{int}}}@f$ and 
%> @f$\delta_{E_{Kn}\in\ensuremath{\mathcal{E}}_{\text{bc}}}@f$ denote the 
%> Kronecker delta, @f$\otimes@f$ denotes the Kronecker product and @f$\mathsf{\Delta}_n@f$
%> is the permutation matrix mapping from the element-based view of the assembly towards the 
%> edge-based view of the hybrid degrees of freedom (see <code>assembleMatEdgePhiIntMuVal()</code>).
%> Which mass matrix is assembled by the function depends on the edges marked
%> by \p markE0T.
%> 
%> The entries of matrix @f$\hat{\mathsf{M}}_\mu \in \mathbb{R}^{\bar{N}\times\bar{N}}@f$ are given by
%> @f[
%> [\hat{\mathsf{M}}_\mu]_{i,j} = \int_0^1 \hat{\mu}_{j}(s) \, \hat{\mu}_{i}(s) \, \text{d}s
%> @f]
%%> @f[
%%> \int_{E_{kn}} \mu_{knj} \, \mu_{kni} \, \text{d}s 
%%> = \ensuremath{|E_{kn}|} \int_0^1 \hat{\mu}_{j} \circ \hat{\beta}_{kn}(s) \, \hat{\mu}_{i} \circ \hat{\beta}_{kn}(s) \, \text{d}s
%%> = \ensuremath{|E_{kn}|} \underbrace{\int_0^1 \hat{\mu}_{j}(s) \, \hat{\mu}_{i}(s) \, \text{d}s}_{\eqqcolon [\hat{\mathsf{M}}_\mu]_{i,j}}\,.
%%> @f]
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
