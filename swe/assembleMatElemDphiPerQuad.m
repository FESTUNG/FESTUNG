% Assembles two matrices, each containing evaluations of a (spatial) derivative
% of one basis function in quadrature points multiplied with the corresponding
% quadrature weight.
%
%===============================================================================
%> @file assembleMatElemDphiPerQuad.m
%>
%> @brief Assembles two matrices, each containing evaluations of a (spatial) 
%>        derivative of one basis function in quadrature points multiplied with
%>        the corresponding quadrature weight.
%===============================================================================
%>
%> @brief Assembles two matrices, each containing evaluations of a (spatial) 
%>        derivative of one basis function in quadrature points multiplied with
%>        the corresponding quadrature weight.
%>
%> The matrices @f$\mathsf{H}^m \in \mathbb{R}^{KN\times KR}@f$ (R is the 
%> number of quadrature points and weights.) are block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{H}^m]_{(k-1)N+i,(k-1)N+j} = 
%>      \partial_{x^m} \varphi_{ki}(q_k^r) w_r  \,.
%> @f]
%> All other entries are zero.
%> For the implementation, the element integrals are backtransformed to the
%> reference triangle @f$\hat{T} = \{(0,0), (1,0), (0,1)\}@f$ using an affine
%> mapping @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$
%> defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,.
%> @f]
%> For the transformation of the gradient holds @f$ \mathbf{\nabla} = 
%> ( \hat{\mathbf{\nabla}} \mathbf{F}_k )^{-T} \hat{\mathbf{\nabla}} @f$,
%> resulting in the component-wise rule in @f$\mathbf{x} \in T_k@f$:
%> @f[
%>   \begin{bmatrix} \partial_{x^1} \\ \partial_{x^2} \end{bmatrix} =
%>   \frac{1}{2|T_k|} \begin{bmatrix}
%>      B_k^{22} \partial_{\hat{x}^1} - B_k^{21} \partial_{\hat{x}^2} \\
%>      B_k^{11} \partial_{\hat{x}^2} - B_k^{12} \partial_{\hat{x}^1}
%>   \end{bmatrix} \,.
%> @f]
%> This allows to assemble the matrices as
%> @f[
%>  \mathsf{H}^1 = B_{k}^{22} [\hat{\mathsf{H}}]_{:,:,1}
%>                   - B_{k}^{21} [\hat{\mathsf{H}}]_{:,:,2} \text{ and }
%>  \mathsf{H}^2 = -B_{k}^{12} [\hat{\mathsf{H}}]_{:,:,1}
%>                   + B_{k}^{11} [\hat{\mathsf{H}}]_{:,:,2} \,,
%> @f]
%> with @f$\hat{\mathsf{H}} \in \mathbb{R}^{N \times R \times 2}@f$ defined as
%> @f[
%>  [\hat{\mathsf{H}}]_{i,r,m} = \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i (\hat{q_r}) \hat(w_r)\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhiPerQuad Local matrix @f$\hat{\mathsf{H}}@f$ as provided
%>                    by <code>integrateRefElemDphiPhiPerQuad()</code>.
%>                    @f$[N \times R \times 2]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatElemDphiPerQuad(g, refElemDphiPerQuad)
K = g.numT;
ret = cell(2,1);
ret{1} = kron(spdiags(g.B(:,2,2), 0, K, K), refElemDphiPerQuad(:,:,1)) - kron(spdiags(g.B(:,2,1), 0, K, K), refElemDphiPerQuad(:,:,2));
ret{2} = kron(spdiags(g.B(:,1,1), 0, K, K), refElemDphiPerQuad(:,:,2)) - kron(spdiags(g.B(:,1,2), 0, K, K), refElemDphiPerQuad(:,:,1));
end % function
