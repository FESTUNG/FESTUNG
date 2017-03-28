% Compute integrals on the reference square, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file integrateRefElemTetraPhiPhi.m
%>
%> @brief Compute integrals on the reference square, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals on the reference square @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions.
%>
%> It computes matrices
%> @f$\hat{\mathsf{M}}^s\in\mathbb{R}^{N\times N}, s\in\{1,2,3\}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{M}}^1]_{i,j} =
%>   \int_{\hat{T}} \hat{\varphi}_i \hat{\varphi}_j \mathrm{d}\hat{\mathbf{x}} \,,
%> @f]
%> @f[
%> [\hat{\mathsf{M}}^2]_{i,j} =
%>   \int_{\hat{T}} \hat{\varphi}_i \hat{\varphi}_j \hat{x}^1 \mathrm{d}\hat{\mathbf{x}} \,,
%> @f]
%> @f[
%> [\hat{\mathsf{M}}^1]_{i,j} =
%>   \int_{\hat{T}} \hat{\varphi}_i \hat{\varphi}_j \hat{x}^2 \mathrm{d}\hat{\mathbf{x}} \,.
%> @f]
%>
%> That way it allows to assemble global matrices on meshes with non-affine
%> transformations.
%>
%> @param  N            The local number of degrees of freedom.
%> @param  qOrd         The order of the quadrature rule.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret  The computed matrices as a cell array @f$[3\times 1]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function hatM = integrateRefElemTetraPhiPhi(N, qOrd, basesOnQuad)
[Q1, Q2, W] = quadRuleTensorProduct(qOrd);
hatM = { zeros(N, N), zeros(N, N), zeros(N, N) };
for i = 1 : N
  for j = 1 : i
    ind = sub2ind([N N], [i j], [j i]);
    phi2Dphi2D = basesOnQuad.phi2D{qOrd}(:, i) .* basesOnQuad.phi2D{qOrd}(:, j);
    hatM{1}(ind) = W * phi2Dphi2D;
    hatM{2}(ind) = W * ( phi2Dphi2D .* Q1' );
    hatM{3}(ind) = W * ( phi2Dphi2D .* Q2' );
  end  % for j
end  % for i
end  % function