% Assembles a matrix, containing integrals of products of two linear 
% Lagrangian basis functions.

%===============================================================================
%> @file assembleMatElemPhiLagrPhiLagr.m
%>
%> @brief Assembles a matrix, containing integrals of products of two  
%>        linear Lagrangian basis functions.
%===============================================================================
%>
%> @brief Assembles a matrix, containing integrals of products of two  
%>        linear Lagrangian basis functions.
%>
%> The matrix @f$\mathsf{M}_L \in \mathbb{R}^{g.numV\times g.numV}@f$ is block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{M}]_{i,j} = \int_{T_k} 
%>      \varphi_{i}^L \varphi_{j}^L \mathrm{d} \mathbf{x} \,.
%> @f]
%> All other entries are zero.
%> For the assembly we define element local matrices 
%> @f$\mathsf{M}_{T_k} \in \mathbb{R}^{N\times N}@f$ as
%> @f[
%>   [\mathsf{M}_{T_k}]_{i,j} = \int_{T_k} 
%>      \varphi_{ki} \varphi_{kj} \mathrm{d} \mathbf{x} \,,
%> @f]
%> and rewrite the global matrix as 
%> @f$\mathsf{M} = \mathrm{diag}(\mathsf{M}_{T_1}, \ldots, \mathsf{M}_{T_K})@f$.
%> In this implementation, the element integrals are backtransformed to the
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
%> This allows to rewrite the local matrices as 
%> @f$\mathsf{M}_{T_k} = 2|T_k| \hat{\mathsf{M}}@f$ with
%> @f[
%>  \hat{\mathsf{M}} = \begin{bmatrix}
%>    \hat{\varphi}_1 \hat{\varphi}_1 & \dots & \hat{\varphi}_1 \hat{\varphi}_N \\
%>    \vdots & \ddots & \vdots \\
%>    \hat{\varphi}_N \hat{\varphi}_1 & \dots & \hat{\varphi}_N \hat{\varphi}_N
%>  \end{bmatrix} \in \mathbb{R}^{N\times N}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemPhiLagrPhiLagr Local matrix @f$\hat{\mathsf{M}}@f$ as provided
%>                    by <code>integrateRefElemPhiLagrPhiLagr()</code>.
%>                    @f$[3 \times 3]@f$
%> @retval ret        The assembled matrix @f$[g.numV \times g.numV]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatElemPhiLagrPhiLagr(g, refElemPhiLagrPhiLagr)
ret = sparse(g.numV, g.numV);
for k = 1 : g.numT
  for i = 1 : 3
    for j = 1 : 3
      ret(g.V0T(k,i), g.V0T(k,j)) = ret(g.V0T(k,i), g.V0T(k,j)) + 2 * g.areaT(k) * refElemPhiLagrPhiLagr(i,j);
    end % for
  end % for
end % for
end % function
