% Assembles a matrix containing integrals of products of a 1D basis 
% function with the (spatial) derivative of a 1D basis function and with a 
% discontinuous 1D coefficient function, divided by the smoothed mesh height.

%===============================================================================
%> @file
%>
%> @brief Assembles a matrix containing integrals of products of a 1D basis 
%>        function with the (spatial) derivative of a 1D basis function and with 
%>        a discontinuous 1D coefficient function, divided by the smoothed mesh 
%>        height.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals of products of a 1D basis 
%>        function with the (spatial) derivative of a 1D basis function and with 
%>        a discontinuous 1D coefficient function, divided by the smoothed mesh 
%>        height.
%>
%> The matrix @f$\bar{\mathsf{G}} \in 
%>   \mathbb{R}^{\overline{K}\overline{N}\times \overline{K}\overline{N}}@f$ 
%> is block diagonal and defined component-wise by
%> @f[
%>   \left[\overline{\mathsf{G}}\right]_{(\overline{k}-1)\overline{N}+i,
%>       (\overline{k}-1)\overline{N}+j} :=
%>     \int_{\overline{T}_{\overline{k}}} \frac{1}{H_s} \partial_{x^1} 
%>     \phi_{\overline{k}i}\, \left( \sum_{l=1}^{\overline{N}} 
%>     \overline{U}_{\overline{k}l} \, \phi_{\overline{k}l} \right)
%>     \, \phi_{\overline{k}j} \, \mathrm{d} x^1\,.
%> @f]
%> All other entries are zero.
%> For the implementation, the element integrals are backtransformed to the
%> reference interval @f$[0,1]@f$ using a mapping 
%> @f$\overline{F}_{\overline{k}}:[0,1]\ni \hat{x}\mapsto x^1\in 
%> \overline{T}_{\overline{k}} @f$.
%>
%> With the chain rule for the partial derivative 
%> @f$\partial_{x^1} \frac{1}{|\overline{T}_{\overline{k}}}@f$ and applying a
%> one-dimensional quadrature rule, this allows to write
%> @f{align*}{
%> \Big[\overline{\mathsf{G}}&\Big]_{(\overline{k}-1)\overline{N}+i,(\overline{k}-1)\overline{N}+j} 
%> := \int_{\overline{T}_{\overline{k}}} \frac{1}{H_s} \partial_{x^1} 
%> \phi_{\overline{k}i}\, \left( \sum_{l=1}^{\overline{N}} \overline{U}_{\overline{k}l} \, 
%> \phi_{\overline{k}l} \right) \, \phi_{\overline{k}j} \, \mathrm{d} x^1
%> = \left|\overline{T}_{\overline{k}}\right| \sum_{l=1}^{\overline{N}} 
%>  \int_0^1 \frac{\overline{U}^1_{\overline{k}{l}} + \overline{U}^2_{\overline{k}{l}}\, \hat{x}}{
%>    H_s \circ \overline{F}_{\overline{k}}(\hat{x})} \,
%> \frac{\partial_{\hat{x}} \hat{\phi}_i(\hat{x})}{\left|\overline{T}_{\overline{k}}\right|}  \,
%>  \hat{\phi}_l(\hat{x}) \, \hat{\phi}_j(\hat{x}) \, \mathrm{d}\hat{x}\\
%> &\approx\sum_{l=1}^{\overline{N}} \sum_{r=1}^R \frac{
%>  \overline{U}^1_{\overline{k}{l}} \; \overbrace{
%>    \omega_r\, \partial_{\hat{x}}\hat{\phi}_i(q_r) \,\hat{\phi}_l(q_r) \,\hat{\phi}_j(q_r)
%>  }^{=: [\hat{\overline{\mathsf{G}}}^1]_{i,j,l,r}} 
%>  + \overline{U}^2_{\overline{k}{l}} \; \overbrace{
%>    \omega_r\, \partial_{\hat{x}}\hat{\phi}_i(q_r) \,\hat{\phi}_l(q_r) \,\hat{\phi}_j(q_r) \, q_r
%>  }^{=: [\hat{\overline{\mathsf{G}}}^2]_{i,j,l,r}}
%>  }{H_s\circ\overline{F}_{\overline{k}}(q_r)}
%>  = \sum_{s=1}^2 \sum_{l=1}^{\overline{N}}
%> \overline{U}^s_{\overline{k}{l}} \sum_{r=1}^{R}
%>   \frac{\left[ \hat{\overline{\mathsf{G}}}^s \right]_{i,j,l,r}}{H_s\circ\overline{F}_{\overline{k}}(q_r)} \,.
%> @f}
%>
%> @param dataDisc    A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$ for scalar or @f$[2 \times 1\text{ cell}/
%>                    [2 \times 2 \text{ cell}@f$ of such matrices for
%>                    vectorial or matrix coefficients, respectively.
%> @param heightQ0T   The smoothed mesh height in each quadrature point of the 1D grid
%>                    @f$[\overline{K} \times R]@f$
%> @param refElemDphiPhiPhiPerQuad Local matrices @f$\hat{\mathsf{G}}^s@f$ as provided
%>                    by <code>integrateRefElem1DDphiPhiPhiPerQuad()</code>.
%>                    @f$[\overline{N} \times \overline{N} \times \overline{N} \times R]@f$
%>
%> @retval ret        The assembled matrix @f$[\overline{K}\overline{N}\times \overline{K}\overline{N}]@f$
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
function ret = assembleMatElem1DDphiPhiFuncDiscHeight(dataDisc, heightQ0T, refElemDphiPhiPhiPerQuad)
[K, N] = size(dataDisc{1}); R = size(heightQ0T, 2);
ret = sparse(K*N, K*N);
invHeightQ0T = 1./heightQ0T;
for s = 1 : 2
  for l = 1 : N
    refElemDphiPhiPhiHeight = zeros(K*N, N);
    for r = 1 : R
      refElemDphiPhiPhiHeight = refElemDphiPhiPhiHeight + kron(invHeightQ0T(:,r), refElemDphiPhiPhiPerQuad{s}(:,:,l,r));
    end % for r
    ret = ret + kronVec(spdiags(dataDisc{s}(:,l), 0, K, K), refElemDphiPhiPhiHeight);
  end % for l
end % for s
end % function