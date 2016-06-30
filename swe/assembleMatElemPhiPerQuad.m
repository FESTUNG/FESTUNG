% Assembles a matrix of one basis function evaluated at all quadrature points 
% of all triangles multiplied with the corresponding quadrature weights.
%
%===============================================================================
%> @file assembleMatElemPhiPerQuad.m
%>
%> @brief Assembles a matrix of one basis function evaluated at all quadrature 
%>				points of all triangles multiplied with the corresponding quadrature 
%>				weights.
%===============================================================================
%>
%> @brief Assembles a matrix of one basis function evaluated at all quadrature 
%>				points of all triangles multiplied with the corresponding quadrature 
%>				weights.
%>
%> The matrix @f$\mathsf{E} \in \mathbb{R}^{KN\times KR}@f$ is block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{E}]_{(k-1)N+i,(k-1)R+r} = \varphi_{ki} (q_r) \, w_r.
%> @f]
%> All other entries are zero.
%> For the assembly we define element local matrices 
%> @f$\mathsf{E}_{T_k} \in \mathbb{R}^{N\times R}@f$ as
%> @f[
%>   [\mathsf{E}_{T_k}]_{i,r} = \varphi_{ki}(q_r) \,,
%> @f]
%> and rewrite the global matrix as 
%> @f$\mathsf{E} = \mathrm{diag}(\mathsf{E}_{T_1}, \ldots, \mathsf{E}_{T_K})@f$.
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
%> @f$\mathsf{E}_{T_k} = 2|T_k| \hat{\mathsf{E}}@f$ with
%> @f[
%>  \hat{\mathsf{E}} = \begin{bmatrix} 
%>    \hat{\varphi}_1 (\hat{q}_1) & \dots & 
%>    \hat{\varphi}_l (\hat{q}_R) \\
%>    \vdots & \ddots & \vdots \\
%>    \hat{\varphi}_N (\hat{q}_1) & 
%>    \dots & \hat{\varphi}_N (\hat{q}_R)
%>  \end{bmatrix} \in \mathbb{R}^{N\times R}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemPhiPerQuad Local tensor @f$\hat{\mathsf{E}}@f$ as provided
%>                    by <code>integrateRefElemPhiPerQuad()</code>.
%>                    @f$[N \times R]@f$
%> @retval ret        The assembled matrix @f$[KN \times KR]@f$
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
function ret = assembleMatElemPhiPerQuad(g, refElemPhiPerQuad)
K = g.numT;
ret = 2 * kron(spdiags(g.areaT, 0, K, K), refElemPhiPerQuad);
end % function
