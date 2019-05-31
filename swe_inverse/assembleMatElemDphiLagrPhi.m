% Assembles two matrices, each containing integrals of products of a basis 
% function with a (spatial) derivative of a linear Lagrange basis function.
%
%===============================================================================
%> @file assembleMatElemDphiLagrPhi.m
%>
%> @brief Assembles two matrices, each containing integrals of products of 
%>				a basis function with a (spatial) derivative of a linear Lagrange
%>				basis function.
%===============================================================================
%>
%> @brief Assembles two matrices, each containing integrals of products of 
%>				a basis function with a (spatial) derivative of a linear Lagrange
%>				basis function.
%>
%> The matrices @f$\mathsf{H}_L^m \in \mathbb{R}^{L\times KN}@f$ are defined 
%> component-wise by
%> @f[
%>   [\mathsf{H}_L^m]_{i,(k-1)N+j} = \int_{T_k} 
%>      \partial_{x^m} \varphi_i^L \varphi_{kj} \mathrm{d} \mathbf{x} \,.
%> @f]
%> where \varphi_i^L is the piecewise linear continuous function whose
%> value in the i-th global grid vertex is one and zero in all others.
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
%>  \mathsf{H}_L^1 = B_{k}^{22} [\hat{\mathsf{H}_L}]_{:,:,1}
%>                   - B_{k}^{21} [\hat{\mathsf{H}_L}]_{:,:,2} \text{ and }
%>  \mathsf{H}_L^2 = -B_{k}^{12} [\hat{\mathsf{H}_L}]_{:,:,1}
%>                   + B_{k}^{11} [\hat{\mathsf{H}_L}]_{:,:,2} \,,
%> @f]
%> with @f$\hat{\mathsf{H}_L} \in \mathbb{R}^{3 \times N \times 2}@f$ defined as
%> @f[
%>  [\hat{\mathsf{H}}]_{i,j,m} = \int_{\hat{T}} \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i^L \hat{\varphi}_j \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiLagrPhi Local matrix @f$\hat{\mathsf{H}}@f$ as provided
%>                    by <code>integrateRefElemDphiLagrPhi()</code>.
%>                    @f$[3 \times N \times 2]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatElemDphiLagrPhi(g, refElemDphiLagrPhi)
K = g.numT; N = size(refElemDphiLagrPhi, 2); L = g.numV;

% Check function arguments that are directly used
validateattributes(refElemDphiLagrPhi, {'numeric'}, {'size', [3 N 2]}, mfilename, 'refElemDphiPhi');

% Assemble matrices
ret = cell(2, 1); ret{1} = sparse(L, K*N); ret{2} = sparse(L, K*N);
for n = 1:3
	markVV0T = sparse(bsxfun(@eq, (1:L)', g.V0T(:,n)'));
	ret{1} = ret{1} + kron(bsxfun(@times, markVV0T, g.B(:,2,2).'), refElemDphiLagrPhi(n,:,1)) ...
									- kron(bsxfun(@times, markVV0T, g.B(:,2,1).'), refElemDphiLagrPhi(n,:,2));
	ret{2} = ret{2} - kron(bsxfun(@times, markVV0T, g.B(:,1,2).'), refElemDphiLagrPhi(n,:,1)) ...
									+ kron(bsxfun(@times, markVV0T, g.B(:,1,1).'), refElemDphiLagrPhi(n,:,2));
end % for
end % function
