% Assembles a matrix containing integrals of products of (spatial)
% derivatives of two basis functions.

%===============================================================================
%> @file ./core/assembleMatElemDphiDPhi.m
%>
%> @brief Assembles a matrix containing integrals of products of (spatial)
%>        derivatives of two basis functions.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals of products of (spatial)
%>        derivatives of two basis functions.
%>
%> The matrix @f$\mathsf{L} \in \mathbb{R}^{KN\times KN}@f$ is block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{L}^m]_{(k-1)N+i,(k-1)N+j} = \int_{T_k} 
%>      \nabla \varphi_{ki} \cdot \nabla \varphi_{kj} \mathrm{d} \mathbf{x} \,.
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
%> This allows to assemble the matrices as described in
%> [FrankKuzminRupp2018].
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemDphiDPhi Local matrix @f$\hat{\mathsf{L}}@f$ as provided
%>                    by <code>integrateRefElemDphiDPhi()</code>.
%>                    @f$[N \times N \times 2]@f$
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Andreas Rupp, 2018.
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
function ret = assembleMatElemDPhiDPhi(g, refElemDPhiDPhi)
K = g.numT;

% triangular mesh with affine-linear mapping
N = size(refElemDPhiDPhi, 1);

% Check function arguments that are directly used
validateattributes(refElemDPhiDPhi, {'numeric'}, {'size', [N N 3]}, mfilename, 'refElemDPhiDPhi');

% Assemble matrices
ret = kron(spdiags( ( g.B(:,1,2) .* g.B(:,1,2) + g.B(:,2,2) .* g.B(:,2,2) ) ./ g.detJ0T, 0, K, K), refElemDPhiDPhi(:,:,1)) ...
      - kron(spdiags( ( g.B(:,1,1) .* g.B(:,1,2) + g.B(:,2,1) .* g.B(:,2,2) ) ./ g.detJ0T, 0, K, K), refElemDPhiDPhi(:,:,2)) ...
      + kron(spdiags( ( g.B(:,1,1) .* g.B(:,1,1) + g.B(:,2,1) .* g.B(:,2,1) ) ./ g.detJ0T, 0, K, K), refElemDPhiDPhi(:,:,3));
end % function
