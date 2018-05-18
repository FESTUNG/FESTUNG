% Assembles a matrix containing products of two one-dimensional basis functions
% with a discrete 1D function, all from the interior of the element, with the 
% normal "vector" and divided by the smoothed mesh height.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing products of two one-dimensional basis 
%>        functions with a discrete 1D function, all from the interior of the
%>        element, with the normal "vector" and divided by the smoothed mesh height.
%===============================================================================
%>
%> @brief Assembles a matrix containing products of two one-dimensional basis 
%>        functions with a discrete 1D function, all from the interior of the
%>        element, with the normal "vector" and divided by the smoothed mesh height.
%>
%> The matrix @f$\overline{\mathsf{P}}_\mathrm{bdr} \in 
%>      \mathbb{R}^{\overline{K}\overline{N}\times \overline{K}\overline{N}}@f$
%> is almost identical to the block diagonal part of
%> swe_2dv/assembleMatV0T1DPhiPhiFuncDiscNuHeight.m 
%> Its entries are given as
%> @f[
%> \left[\overline{\mathsf{P}}_\mathrm{bdr}\right]_{(\overline{k}-1)\overline{N}+i,
%>     (\overline{k}-1)\overline{N}+j} :=
%>  \sum_{a^1_{\overline{k}\overline{n}}\in
%>     \partial \overline{T}_{\overline{k}}\cap \overline{\mathcal{V}}_\Omega}
%>  \frac{\nu_{\overline{k}\overline{n}}^1}{H_s\left(a^1_{\overline{k}\overline{n}}\right)}
%>  \sum_{l=1}^{\overline{N}} \overline{U}_{\overline{k}l}\left(a^1_{\overline{k}\overline{n}}\right)\;
%>  \phi_{\overline{k}i} \left(a^1_{\overline{k}\overline{n}}\right) \, 
%>  \phi_{\overline{k}l}\left(a^1_{\overline{k}\overline{n}}\right) \, 
%>  \phi_{\overline{k}j}\left(a^1_{\overline{k}\overline{n}}\right) 
%> =: \frac{1}{2} \sum_{a^1_{\overline{k}\overline{n}}\in
%>     \partial \overline{T}_{\overline{k}}\cap \overline{\mathcal{V}}_\Omega}
%>  \frac{\nu_{\overline{k}\overline{n}}^1}{H_s\left(a^1_{\overline{k}\overline{n}}\right)}
%>  \sum_{s=1}^2 \sum_{l=1}^{\overline{N}} \overline{U}^s_{\overline{k}l} 
%>  \left[\hat{\bar{\mathsf{P}}}^{s,\text{diag}}\right]_{i,j,l,\overline{n}} \,,
%> @f]
%> with @f$\nu_{\overline{k}^-\overline{n}^-}@f$ the "normal" (i.e., @f$\pm 1@f$)
%> and @f$\overline{U}_{\overline{k}l}@f$ the @f$(k,l)@f$-th entry of the 
%> representation matrix of the discrete function @f$\overline{u}_\Delta@f$.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a 1D triangulation (see 
%>                    <code>generateGridData1D()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param dataDisc    A representation matrix of the discrete function 
%>                    @f$\overline{u}_\Delta(x^1)@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc1D()</code>
%>                    @f$[\overline{K} \times \overline{N}]@f$
%> @param heightV0T   The smoothed mesh height in each node of the 1D grid
%>                    @f$[\overline{K} \times 2]@f$
%> @param  markV0T    <code>logical</code> arrays that mark each elements
%>                    vertices on which the matrix blocks should be
%>                    assembled @f$[\overline{K} \times 2]@f$
%> @param refEdgePhiIntPhiIntPhiInt  Local matrices 
%>                    @f$\hat{\bar{\mathsf{P}}}^\text{diag}@f$ as provided
%>                    by <code>computePhiIntPhiIntPhiIntV0T1D()</code>.
%>                    @f$[2 \times 1 \text{ cell}]@f$
%>
%> @retval ret        The assembled matrix @f$[\bar{K}\bar{N} \times \bar{K}\bar{N}]@f$
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
function ret = assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(g, dataDisc, heightV0T, markV0T, refEdgePhiIntPhiIntPhiInt)
[K,N] = size(dataDisc{1});
ret = sparse(K*N, K*N);
for n = 1 : 2
  nuV0THeight = markV0T(:,n) .* g.nuV0T(:,n) ./ heightV0T(:,n);
  for s = 1 : 2
    for l = 1 : N
      ret = ret + kron(spdiags(nuV0THeight .* dataDisc{s}(:,l), 0, K, K), ...
                       refEdgePhiIntPhiIntPhiInt{s}(:,:,l,n));
    end % for l
  end % for s
end % for n