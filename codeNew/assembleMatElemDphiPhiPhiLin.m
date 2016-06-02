% Assembles two matrices, each containing integrals of products of a basis 
% function and a linear basis function with a (spatial) derivative of a basis 
% function.
%
%===============================================================================
%> @file assembleMatElemDphiPhiPhiLin.m
%>
%> @brief Assembles two matrices, each containing integrals of products of a 
%>        basis function and a linear basis function with a (spatial) derivative
%>				of a basis function.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{H}^m, m \in \{1,2\}@f$
%>        containing integrals of products of a basis function and a linear 
%>				basis function with a (spatial) derivative of a basis function.
%>
%> The matrices @f$\mathsf{H}^m \in \mathbb{R}^{KN\times KN}@f$ are block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{H}^m]_{(k-1)N+i,(k-1)N+j} = sum_{l=1}^3 zbDisc(k,l) \int_{T_k} 
%>      \partial_{x^m} \varphi_{ki} \varphi_{kj} 
%>			\varphi_{kl} \mathrm{d} \mathbf{x} \,.
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
%> @f[ sum_{l=1}^3
%>  \mathsf{H}^1 = zbDisc(k,l) \left( B_{k}^{22} [\hat{\mathsf{H}}]_{:,:,l,1}
%>                   - B_{k}^{21} [\hat{\mathsf{H}}]_{:,:,l,2} \right) \text{ and }
%>  \mathsf{H}^2 = zbDisc(k,l) \left( -B_{k}^{12} [\hat{\mathsf{H}}]_{:,:,l,1}
%>                   + B_{k}^{11} [\hat{\mathsf{H}}]_{:,:,l,2} \right) \,,
%> @f]
%> with @f$\hat{\mathsf{H}} \in \mathbb{R}^{N \times N \times 3\times 2}@f$ defined as
%> @f[
%>  [\hat{\mathsf{H}}]_{i,j,l,m} = \int_{\hat{T}} \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i \hat{\varphi}_j \hat{\varphi}_l \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhiPhiLin Local matrix @f$\hat{\mathsf{H}}@f$ as provided
%>                    by <code>integrateRefElemDphiPhiPhiLin()</code>.
%>                    @f$[N \times N \times 3\times 2]@f$
%> @param dataDiscLin	A representation of the discrete function ,e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
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
function ret = assembleMatElemDphiPhiPhiLin(g, refElemDphiPhiPhiLin, dataDiscLin)
K = g.numT; N = size(refElemDphiPhiPhiLin, 1);

% Check function arguments that are directly used
validateattributes(refElemDphiPhiPhiLin, {'numeric'}, {'size', [N N 3 2]}, mfilename, 'refElemDphiPhi');

% Assemble matrices
ret = cell(2, 1); ret{1} = sparse(K*N, K*N); ret{2} = sparse(K*N, K*N);
for l = 1 : 3
	ret{1} = ret{1} + kron(spdiags(g.B(:,2,2) .* dataDiscLin(:,l), 0,K,K), refElemDphiPhiPhiLin(:,:,l,1)) ...
									- kron(spdiags(g.B(:,2,1) .* dataDiscLin(:,l), 0,K,K), refElemDphiPhiPhiLin(:,:,l,2));
	ret{2} = ret{2} - kron(spdiags(g.B(:,1,2) .* dataDiscLin(:,l), 0,K,K), refElemDphiPhiPhiLin(:,:,l,1)) ...
									+ kron(spdiags(g.B(:,1,1) .* dataDiscLin(:,l), 0,K,K), refElemDphiPhiPhiLin(:,:,l,2));
end % for
end % function
