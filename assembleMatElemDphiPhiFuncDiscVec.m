% Assembles two matrices, each containing integrals of products of a basis 
% function with a (spatial) derivative of a basis function and with a 
% discontinuous vector-wise coefficient function.
%
%===============================================================================
%> @file
%>
%> @brief Assembles two matrices, each containing integrals of products of a 
%>        basis function with a (spatial) derivative of a basis function
%>        and with a discontinuous vector-wise coefficient function.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{G}^m, m \in \{1,2\}@f$
%>        containing integrals of products of a basis function with a (spatial)
%>        derivative of a basis function and with a discontinuous vector-wise
%>        coefficient function.
%>
%> The matrices @f$\mathsf{G}^m \in \mathbb{R}^{KN\times KN}@f$ are block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{G}^m]_{(k-1)N+i,(k-1)N+j} = \sum_{l=1}^N U^m_{kl}(t) \int_{T_k} 
%>     \partial_{x^m}\varphi_{ki}\varphi_{kl}\varphi_{kj}\mathrm{d}\mathbf{x}\,.
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
%> This allows to write
%> @f[
%>   \int_{T_k} \partial_{x^1} \varphi_{ki} \varphi_{kl} \varphi_{kj} =
%>   B_k^{22}[\hat{\mathsf{G}}]_{i,j,l,1}-B_k^{21}[\hat{\mathsf{G}}]_{i,j,l,2}
%>   \text{ and }
%>   \int_{T_k} \partial_{x^2} \varphi_{ki} \varphi_{kl} \varphi_{kj} =
%>   -B_k^{12}[\hat{\mathsf{G}}]_{i,j,l,1}+B_k^{11}[\hat{\mathsf{G}}]_{i,j,l,2}\,,
%> @f]
%> with @f$\hat{\mathsf{G}} \in \mathbb{R}^{N\times N\times N\times 2}@f$
%> given as 
%> @f[
%>  [\hat{\mathsf{G}}]_{i,j,l,m} = \int_{\hat{T}} 
%>    \partial_{\hat{x}^m} \hat{\varphi}_i \hat{\varphi}_l \hat{\varphi}_j \,.
%> @f]
%> Now we can build local matrices
%> @f$\mathsf{G}_{T_k}^m \in \mathbb{R}^{N\times N}@f$ as
%> @f[
%>  \mathsf{G}_{T_k}^1 = \sum_{l=1}^N U^1_{kl}(t) \left(
%>    B_k^{22}[\hat{\mathsf{G}}]_{:,:,l,1} 
%>    - B_k^{21}[\hat{\mathsf{G}}]_{:,:,l,2} \right)
%> @f]
%> and
%> @f[
%>  \mathsf{G}_{T_k}^2 = \sum_{l=1}^N U^2_{kl}(t) \left(
%>    -B_k^{12}[\hat{\mathsf{G}}]_{:,:,l,1} 
%>    + B_k^{11}[\hat{\mathsf{G}}]_{:,:,l,2} \right) \,.
%> @f]
%> With @f$\mathsf{G}^m = 
%> \mathrm{diag}(\mathsf{G}_{T_1}^m,\ldots,\mathsf{G}_{T_K}^m)@f$
%> we can vectorize the assembly using the Kronecker product.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhiPhi Local matrix @f$\hat{\mathsf{G}}@f$ as provided
%>                    by <code>integrateRefElemDphiPhiPhi()</code>.
%>                    @f$[N \times N \times N \times 2]@f$
%> @param dataDisc1   A representation of the first component of the
%>                    discrete function, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code> @f$[K \times N]@f$
%> @param dataDisc2   A representation of the second component of the
%>                    discrete function. @f$[K \times N]@f$
%> @param markElem    (optional) <code>logical</code> arrays to provide the 
%>                    elements of the grid for which the computation is done.
%>                    @f$[K \times 1]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2016.
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
function ret = assembleMatElemDphiPhiFuncDiscVec(g, refElemDphiPhiPhi, dataDisc1, dataDisc2, markElem)
[K, dataN] = size(dataDisc1);
N = size(refElemDphiPhiPhi,1);

if nargin == 4
  markElem = true(K, 1);
end % if

% Check function arguments that are directly used
validateattributes(dataDisc1, {'numeric'}, {'size', [g.numT dataN]}, mfilename, 'dataDisc1');
validateattributes(dataDisc2, {'numeric'}, {'size', [g.numT dataN]}, mfilename, 'dataDisc2');
validateattributes(refElemDphiPhiPhi, {'numeric'}, {'size', [N N dataN 2]}, mfilename, 'refElemDphiPhiPhi');
validateattributes(markElem, {'logical'}, {'size', [K 1]}, mfilename, 'markElem');

Ke = sum(markElem);

% Assemble matrices
ret = cell(2, 1);  ret{1} = sparse(Ke*N, Ke*N);  ret{2} = sparse(Ke*N, Ke*N);
for l = 1 : dataN
  ret{1} = ret{1} + kron(spdiags(dataDisc1(markElem,l) .* g.B(markElem,2,2), 0,Ke,Ke), refElemDphiPhiPhi(:,:,l,1)) ...
                  - kron(spdiags(dataDisc1(markElem,l) .* g.B(markElem,2,1), 0,Ke,Ke), refElemDphiPhiPhi(:,:,l,2));
  ret{2} = ret{2} - kron(spdiags(dataDisc2(markElem,l) .* g.B(markElem,1,2), 0,Ke,Ke), refElemDphiPhiPhi(:,:,l,1)) ...
                  + kron(spdiags(dataDisc2(markElem,l) .* g.B(markElem,1,1), 0,Ke,Ke), refElemDphiPhiPhi(:,:,l,2));
end % for
end % function
