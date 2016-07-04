% Assembles a matrix, containing integrals of products of two basis functions
% with a discontinuous coefficient function.

%===============================================================================
%> @file assembleMatElemPhiPhiFuncDisc.m
%>
%> @brief Assembles a matrix, containing integrals of products of two basis 
%>        functions with a discontinuous coefficient function.
%===============================================================================
%>
%> @brief Assembles a matrix, containing integrals of products of two basis 
%>        functions with a discontinuous coefficient function.
%>
%> The matrix @f$\mathsf{D} \in \mathbb{R}^{KN\times KN}@f$ is block
%> diagonal and defined component-wise by
%> @f[
%>   [\mathsf{D}]_{(k-1)N+i,(k-1)N+j} = \sum_{l=1}^{N_\mathrm{data}} {fc}_{kl}
%>      \int_{T_k} \varphi_{ki}\varphi_{kl}\varphi_{kj}\mathrm{d}\mathbf{x}\,.
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
%> For the assembly we define element local matrices 
%> @f$\mathsf{D}_{T_k} \in \mathbb{R}^{N\times N \times 3}@f$ as
%> @f[
%>   [\mathsf{D}_{T_k}]_{i,j,l} = \int_{T_k} 
%>      \varphi_{ki} \varphi_{kj} \varphi_{kl} \mathrm{d} \mathbf{x} \,,
%> @f]
%> and rewrite the global matrix as 
%> @f$\mathsf{D} = \mathrm{diag}(\mathsf{D}_{T_1}, \ldots, \mathsf{D}_{T_K})@f$.
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
%> @f$\mathsf{D}_{T_k} = 2|T_k| \hat{\mathsf{D}}@f$ with
%> @f[
%>  \hat{\mathsf{D}} = \begin{bmatrix} \sum_{l=1}^{N_\mathrm{data}} 
%>    \hat{\varphi}_l \hat{\varphi}_1 \hat{\varphi}_1 & \dots & 
%>    \hat{\varphi}_l \hat{\varphi}_1 \hat{\varphi}_N \\
%>    \vdots & \ddots & \vdots \\
%>    \hat{\varphi}_l \hat{\varphi}_N \hat{\varphi}_1 & 
%>    \dots & \hat{\varphi}_l \hat{\varphi}_N \hat{\varphi}_N
%>  \end{bmatrix} \in \mathbb{R}^{N\times N}\,.
%> @f]
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param refElemPhiPhiPhi Local matrix @f$\hat{\mathsf{D}}@f$ as provided
%>                    by <code>integrateRefElemPhiPhiPhi()</code>.
%>                    @f$[N \times N \times {N_\mathrm{data}}]@f$
%> @param dataDisc    A representation of the discrete function ,e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times {N_\mathrm{data}}]@f$
%> @retval ret        The assembled matrix @f$[KN \times KN]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatElemPhiPhiFuncDisc(g, refElemPhiPhiPhi, dataDisc)
[K, dataN] = size(dataDisc); 
[N1, N2, ~] = size(refElemPhiPhiPhi);

% Check function arguments that are directly used
validateattributes(dataDisc, {'numeric'}, {'size', [g.numT dataN]}, mfilename, 'dataDisc');
validateattributes(refElemPhiPhiPhi, {'numeric'}, {'size', [N1 N2 dataN]}, mfilename, 'refElemPhiPhiPhi');

ret = sparse(K*N1,K*N2);
for l = 1 : dataN
  ret = ret + 2 * kron(spdiags(g.areaT .* dataDisc(:,l), 0, K, K), refElemPhiPhiPhi(:,:,l));
end % for
end % function
