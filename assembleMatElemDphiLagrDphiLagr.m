% Assembles a matrix containing integrals of products of (spatial) derivatives 
% of two linear Lagrange basis function.
%
%===============================================================================
%> @file assembleMatElemDphiLagrDphiLagr.m
%>
%> @brief Assembles a matrix containing integrals of products of (spatial) 
%>        derivatives of two linear Lagrange basis function.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals of products of (spatial) 
%>        derivatives of two linear Lagrange basis function.
%>
%> The matrix @f$\mathsf{A}_L \in \mathbb{R}^{L\times L}@f$ is defined 
%> component-wise by
%> @f[
%>   [\mathsf{A}_L]_{i,j} = \sum_{T_k \in \Tau_h}\int_{T_k} 
%>      \partial_{x^1} \varphi_i^L \partial_{x^1}\varphi_j^L + 
%>      \partial_{x^2} \varphi_i^L \partial_{x^2}\varphi_j^L \mathrm{d} \mathbf{x} \,.
%> @f]
%> where \varphi_i^L is the piecewise linear continuous function whose
%> value in the i-th global grid vertex is one and zero in all others.
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
%>
%> For the assembly we refer to the book by Knabner and Angerman:
%> Numerical Methods for Elliptic and Parabolic Partial Differential Equations
%>
%> and use @f$\hat{\mathsf{A}_L} \in \mathbb{R}^{3 \times 3 \times 2 \times 2}@f$ defined as
%> @f[
%>  [\hat{\mathsf{A}}]_{i,j,m_1,m_2} = \int_{\hat{T}} \partial_{\hat{x}^{m_1}} 
%>    \hat{\varphi}_i^L \partial_{\hat{x}^{m_2}}\hat{\varphi}_j \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiLagrDphiLagr Local matrix @f$\hat{\mathsf{A}}@f$ as provided
%>                    by <code>integrateRefElemDphiLagrDphiLagr()</code>.
%>                    @f$[3 \times 3 \times 2 \times 2]@f$
%> @retval ret        The assembled matrix @f$[L, L]@f$
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
function ret = assembleMatElemDphiLagrDphiLagr(g, refElemDphiLagrDphiLagr)
K = g.numT; L = g.numV;

% Check function arguments that are directly used
validateattributes(refElemDphiLagrDphiLagr, {'numeric'}, {'size', [3 3 2 2]}, mfilename, 'refElemDphiLagrDphiLagr');

% Assemble matrices
ret = sparse(L, L);
for k = 1 : K
  for i = 1 : 3
    for j = 1 : 3
      ret(g.V0T(k,i),g.V0T(k,j)) = ret(g.V0T(k,i),g.V0T(k,j)) + 0.5 / g.areaT(k) * ( g.B(k,:,2)*g.B(k,:,2)' * refElemDphiLagrDphiLagr(i,j,1,1) ...
                                                                                   - g.B(k,:,1)*g.B(k,:,2)' * (refElemDphiLagrDphiLagr(i,j,1,2)+refElemDphiLagrDphiLagr(i,j,2,1)) ...
                                                                                   + g.B(k,:,1)*g.B(k,:,1)' * refElemDphiLagrDphiLagr(i,j,2,2) );
    end % for
  end % for
end % for
end % function
