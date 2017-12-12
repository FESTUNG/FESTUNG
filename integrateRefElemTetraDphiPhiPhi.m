% Compute integrals on the reference square, whose integrands consist of all
% permutations of two basis functions with one of the (spatial) derivatives
%  of a basis function.

%===============================================================================
%> @file integrateRefElemTetraDphiPhiPhi.m
%>
%> @brief Compute integrals on the reference square, whose integrands consist 
%>        of all permutations of two basis functions with one of the (spatial)
%>        derivatives of a basis function.
%===============================================================================
%>
%> @brief Compute integrals on the reference square @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions
%>        with one of the (spatial) derivatives of a basis function.
%>
%> It computes multidimensional arrays
%> @f$\hat{\mathsf{G}}^s \in \mathbb{R}^{N \times N \times N \times 2}, s\in\{1,2,3\}@f$
%> defined by
%> @f[
%>  [\hat{\mathsf{G}}^1]_{i,j,l,m} = \int_{\hat{T}} \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i \hat{\varphi}_j \hat{\varphi}_l \mathrm{d} \hat{\mathbf{x}}\,,
%> @f]
%> @f[
%>  [\hat{\mathsf{G}}^2]_{i,j,m} = \int_{\hat{T}} \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i \hat{\varphi}_j \hat{\varphi}_l \hat{x}^1 \mathrm{d} \hat{\mathbf{x}}\,,
%> @f]
%> @f[
%>  [\hat{\mathsf{G}}^3]_{i,j,m} = \int_{\hat{T}} \partial_{\hat{x}^m} 
%>    \hat{\varphi}_i \hat{\varphi}_j \hat{\varphi}_l \hat{x}^2 \mathrm{d} \hat{\mathbf{x}}\,,
%> @f]
%>
%> That way it allows to assemble global matrices on meshes with non-affine
%> transformations.
%>
%> @param  N            The local number of degrees of freedom.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D and gradPhi2D.
%> @param  qOrd         The order of the quadrature rule.
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
function ret = integrateRefElemTetraDphiPhiPhi(N, basesOnQuad, qOrd)
[Q1, Q2, W] = quadRuleTensorProduct(qOrd);
ret = { zeros(N, N, N, 2), zeros(N, N, N, 2), zeros(N, N, N, 2) };
for m = 1 : 2
  for j = 1 : N
    for l = 1 : j
      for i = 1 : N
        ind = sub2ind([N, N, N, 2], [i i], [j l], [l j], [m m]);
        gradPhi2Dphi2Dphi2D = basesOnQuad.gradPhi2D{qOrd}(:,i,m) .* basesOnQuad.phi2D{qOrd}(:,j) .* basesOnQuad.phi2D{qOrd}(:,l);
        ret{1}(ind) = W * gradPhi2Dphi2Dphi2D;
        ret{2}(ind) = W * (gradPhi2Dphi2Dphi2D .* Q1');
        ret{3}(ind) = W * (gradPhi2Dphi2Dphi2D .* Q2');
      end  % for i
    end  % for l
  end  % for j
end  % for m
end  % function