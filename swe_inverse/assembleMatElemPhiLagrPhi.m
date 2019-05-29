% Assembles a matrix containing integrals of products of a linear Lagrange basis 
% function with a DG basis function.
%
%===============================================================================
%> @file assembleMatElemPhiLagrPhi.m
%>
%> @brief Assembles a matrix containing integrals of products of a linear 
%>        Lagrange basis function with a DG basis function.
%===============================================================================
%>
%> @brief Assembles a matrix containing integrals of products of a linear 
%>        Lagrange basis function with a DG basis function.
%>
%> The matrix @f$\mathsf{M}_m \in \mathbb{R}^{g.numV\times KN}@f$ is defined 
%> component-wise by
%> @f[
%>   [\mathsf{M}_m]_{i,(k-1)N+j} = \int_{T_k} 
%>      \varphi_i^L \varphi_{kj} \mathrm{d} \mathbf{x} \,.
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
%> This allows to assemble the matrix as
%> @f[
%>  \mathsf{M}_m = 2 |T_k| [\hat{\mathsf{M}_m}]\,,
%> @f]
%> with @f$\hat{\mathsf{M}_m} \in \mathbb{R}^{3 \times N}@f$ defined as
%> @f[
%>  [\hat{\mathsf{M}_m}]_{i,j} = \int_{\hat{T}} 
%>    \hat{\varphi}_i^L \hat{\varphi}_j \mathrm{d} \hat{\mathbf{x}}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemPhiLagrPhi Local matrix @f$\hat{\mathsf{H}}@f$ as provided
%>                    by <code>integrateRefElemPhiLagrPhi()</code>.
%>                    @f$[3 \times N]@f$
%> @retval ret        The assembled matrix @f$[g.numV \times KN]@f$
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
function ret = assembleMatElemPhiLagrPhi(g, refElemPhiLagrPhi)
K = g.numT; N = size(refElemPhiLagrPhi, 2); L = g.numV;

% Check function arguments that are directly used
validateattributes(refElemPhiLagrPhi, {'numeric'}, {'size', [3 N]}, mfilename, 'refElemDphiPhi');

% Assemble matrices
ret = sparse(L, K*N);
for i = 1:3
	markVV0T = sparse(bsxfun(@eq, (1:L)', g.V0T(:,i)'));
	ret = ret + kron(bsxfun(@times, markVV0T, 2 * g.areaT'), refElemPhiLagrPhi(i,:));
end % for
end % function
