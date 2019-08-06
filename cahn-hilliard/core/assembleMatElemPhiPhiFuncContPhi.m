% Assembles a block diagonal matrix containing integrals of products of
% basis functions with a non-linear function depending on the solution

%===============================================================================
%> @file ./core/assembleMatElemPhiPhiFuncContPhi.m
%>
%> @brief Assembles a block diagonal matrix containing integrals of products
%>        of basis functions with a non-linear function depending on the solution.
%===============================================================================
%>
%> @brief Assembles a block diagonal matrix containing integrals of products
%>        of basis functions with a non-linear function depending on the solution.
%>
%> The matrix @f$\mathsf{M}_{\Psi''_+} \in \mathbb{R}^{KN\times KN}@f$ is
%> block diagonal and defined component-wise by
%> @f[
%>   [\mathsf{M}_{\Psi''_+}]_{(k-1)N+i,(k-1)N+j} = \int_{T_k} 
%>     \varphi_{ki}\,\varphi_{kj}\; \Psi''_+\left(\sum_{s=1}^N c_{ks} \varphi_{ks}\right)
%>     \mathrm{d}\mathbf{x}\,.
%> @f]
%> All other entries are zero.
%> For the implementation, the element integrals are backtransformed to the
%> reference element @f$\hat{T}@f$ using a mapping 
%> @f$\mathbf{F}_k:\hat{T}\ni\hat{\mathbf{x}}\mapsto\mathbf{x}\in T_k@f$
%> defined as
%> @f[
%> \mathbf{F}_k (\hat{\mathbf{x}}) = 
%>   \mathsf{{B}}_k \hat{\mathbf{x}} + \hat{\mathbf{a}}_{k1}
%>   \text{ with }
%> \mathbb{R}^{2\times2} \ni \mathsf{{B}}_k =
%>   \left[ \hat{\mathbf{a}}_{k2} - \hat{\mathbf{a}}_{k1} | 
%>          \hat{\mathbf{a}}_{k3} - \hat{\mathbf{a}}_{k1} \right] \,.
%> @f]
%>
%> This allows to assemble the matrices as described in
%> [FrankKuzminRupp2018].
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see, e.g.,
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @param dataDisc    A representation of the discrete function, e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$ for scalar or @f$[2 \times 1\text{ cell}/
%>                    [2 \times 2 \text{ cell}@f$ of such matrices for
%>                    vectorial or matrix coefficients, respectively.
%> @param func        A function handle for the continuous function.
%> @retval ret        The assembled matrices @f$[2 \times 1 \text{ cell}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = assembleMatElemPhiPhiFuncContPhi(g, basesOnQuad, dataDisc, func, hatQuadPhiPhi)
% Extract relevant grid information
[K, N] = size(dataDisc);
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);
[~,~,W] = quadRule2D(qOrd);
R = length(W);

% Compute function values on quadrature points
funcVal = zeros(K,size(basesOnQuad.phi2D{qOrd}(:, 1),1));
for l = 1 : N
  funcVal = funcVal + kron(dataDisc(:,l) , basesOnQuad.phi2D{qOrd}(:, l).');
end
funcVal = func(funcVal);

if nargin < 5
  % Compute value of phi_i * phi_j on quadrature point r
  hatQuadPhiPhi = zeros(size(basesOnQuad.phi2D{qOrd}, 2), size(basesOnQuad.phi2D{qOrd}, 2), size(basesOnQuad.phi2D{qOrd}, 1));
  for r = 1 : R
    hatQuadPhiPhi(:,:,r) = W(r) * kron( basesOnQuad.phi2D{qOrd}(r,:), basesOnQuad.phi2D{qOrd}(r,:).' );
  end % for r = 1 : R
end % if nargin < 5

% Integrate phi_i * phi_j * function_value with respect to elements
ret = zeros(K*N,N);
for r = 1 : R
  ret = ret + kron( g.detJ0T .* funcVal(:,r), hatQuadPhiPhi(:,:,r) );
end
ret = kronVec(speye(K,K), ret);

end % function assembleMatElemPhiPhiFuncContPhi