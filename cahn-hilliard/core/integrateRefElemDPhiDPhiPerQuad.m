% AUTOR: TIM ROITH
% Assembles two matrices, each containing integrals of products of a basis 
% function with a (spatial) derivative of a basis function and with a 
% discontinuous coefficient function.

%===============================================================================
%> @file assembleMatElemDphiPhiFuncDisc.m
%>
%> @brief Assembles two matrices, each containing integrals of products of a 
%>        basis function with a (spatial) derivative of a basis function
%>        and with a discontinuous coefficient function.
%===============================================================================
%>
%> @brief Assembles matrices @f$\mathsf{G}^m, m \in \{1,2\}@f$
%>        containing integrals of products of a basis function with a (spatial)
%>        derivative of a basis function and with a discontinuous coefficient 
%>        function.
%>
%> The matrices @f$\mathsf{G}^m \in \mathbb{R}^{KN\times KN}@f$ are block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{G}^m]_{(k-1)N+i,(k-1)N+j} = \sum_{l=1}^N D_{kl}(t) \int_{T_k} 
%>     \partial_{x^m}\varphi_{ki}\varphi_{kl}\varphi_{kj}\mathrm{d}\mathbf{x}\,.
%> @f]
%> All other entries are zero.
%> For the implementation, the element integrals are backtransformed to the
%> reference element @f$\hat{T}@f$ using a mapping 
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$.
%>
%> In the case of an affine mapping, it is defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,.
%> @f]
%>
%> For the transformation of the gradient holds @f$ \mathbf{\nabla} = 
%> ( \hat{\mathbf{\nabla}} \mathbf{F}_k )^{-T} \hat{\mathbf{\nabla}} =
%> (\mathbf{J}_k(\vec{x}))^{-T} \hat{\mathbf{\nabla}} @f$,
%> with @f$\mathbf{J}_k(\vec{x})@f$ the Jacobian of the mapping.
%>
%> This allows to write
%> @f[
%>   \int_{T_k} \partial_{x^m} \varphi_{ki} \varphi_{kl} \varphi_{kj} =
%>   \sum_{s=1}^S [\mathbf{J}_k^s]_{3-m,3-m} [\hat{\mathsf{G}}^s]_{i,j,l,m}
%>    - [\mathbf{J}_k^s]_{3-m,m} [\hat{\mathsf{G}}^s]_{i,j,l,3-m}\,,
%> @f]
%> with @f$S\in\mathbb{N}@f$ the number of parts in the Jacobian. 
%> For affine mappings, the Jacobian is constant for each element, thus @f$S = 1@f$.
%> For trapezoidals, we have @f$S = 3@f$ corresponding to a constant, an
%> @f$x^1@f$- and an @f$x^2@f$-dependent part.
%> @f$\hat{\mathsf{G}} \in \mathbb{R}^{N\times N\times N\times 2}@f$
%> is given as 
%> @f[
%>  [\hat{\mathsf{G}}]_{i,j,l,m} = \int_{\hat{T}} 
%>    \partial_{\hat{x}^m} \hat{\varphi}_i \hat{\varphi}_l \hat{\varphi}_j \,.
%> @f]
%>
%> Now we can build local matrices
%> @f$\mathsf{G}_{T_k}^m \in \mathbb{R}^{N\times N}@f$ as
%> @f[
%>  \mathsf{G}_{T_k}^m = \sum_{l=1}^N D_{kl}(t) \sum_{s=1}^S \left(
%>     [\mathbf{J}_k^s]_{3-m,3-m} [\hat{\mathsf{G}}^s]_{i,j,l,m}
%>    - [\mathbf{J}_k^s]_{3-m,m} [\hat{\mathsf{G}}^s]_{i,j,l,3-m} \right)\,.
%> @f]
%> With @f$\mathsf{G}^m = 
%> \mathrm{diag}(\mathsf{G}_{T_1}^m,\ldots,\mathsf{G}_{T_K}^m)@f$
%> we can vectorize the assembly using the Kronecker product.
%>
%> Optionally, the discrete function given by coefficients @f$D_{kl}@f$ can
%> be different for each spatial dimension @f$m@f$ or even a matrix, thus
%> introducing an additional sum for each @f$m@f$ over the components of
%> the matrix-vector-product.
%>
%> The approximation order of the coefficient function can be different
%> from the local degrees of freedom.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g.,
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiPhiPhi Local matrix @f$\hat{\mathsf{G}}@f$ as provided
%>                    by <code>integrateRefElemDphiPhiPhi()</code> or
%>                    <code>integrateRefElemTetraDphiPhiPhi()</code>.
%>                    @f$[N \times N \times N \times 2]@f$
%> @param dataDisc    A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$ for scalar or @f$[2 \times 1\text{ cell}/
%>                    [2 \times 2 \text{ cell}@f$ of such matrices for
%>                    vectorial or matrix coefficients, respectively.
%> @param markElem    (optional) <code>logical</code> arrays to provide the 
%>                    elements of the grid for which the computation is done.
%>                    @f$[K \times 1]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Hennes Hajduk, 2016.
%> @author Balthasar Reuter, 2017.
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
function hatQuadDPhiDPhi = integrateRefElemDPhiDPhiPerQuad(N, basesOnQuad)
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);
[~,~,W] = quadRule2D(qOrd);
R = length(W);
%squeeze
hatQuadDPhiDPhi = zeros(size(basesOnQuad.gradPhi2D{qOrd}, 2), size(basesOnQuad.gradPhi2D{qOrd}, 2),3, R);
for r = 1 : R
    hatQuadDPhiDPhi(:,:,1,r) = W(r) * basesOnQuad.gradPhi2D{qOrd}(r,:,1).' * basesOnQuad.gradPhi2D{qOrd}(r,:,1);
    hatQuadDPhiDPhi(:,:,2,r) = W(r) * (basesOnQuad.gradPhi2D{qOrd}(r,:,1).' * basesOnQuad.gradPhi2D{qOrd}(r,:,2) +...
        basesOnQuad.gradPhi2D{qOrd}(r,:,2).' * basesOnQuad.gradPhi2D{qOrd}(r,:,1));
    hatQuadDPhiDPhi(:,:,3,r) = W(r) * basesOnQuad.gradPhi2D{qOrd}(r,:,2).' * basesOnQuad.gradPhi2D{qOrd}(r,:,2);
end % for r = 1 : R
end



