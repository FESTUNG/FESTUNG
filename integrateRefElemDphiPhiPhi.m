% Compute integrals on the reference triangle, whose integrands consist of all
% permutations of two basis functions with one of the (spatial) derivatives of a
% basis function.

%===============================================================================
%> @file integrateRefElemDphiPhiPhi.m
%>
%> @brief Compute integrals on the reference triangle, whose integrands consist 
%>        of all permutations of two basis functions with one of the (spatial)
%>        derivatives of a basis function.
%===============================================================================
%>
%> @brief Compute integrals on the reference triangle @f$\hat{T}@f$,
%>        whose integrands consist of all permutations of two basis functions
%>        with one of the (spatial) derivatives of a basis function.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{G}} \in \mathbb{R}^{N\times N\times N\times 2}@f$
%> defined by
%> @f[
%>  [\hat{\mathsf{G}}]_{i,j,l,m} = \int_{\hat{T}} 
%>    \partial_{\hat{x}^m} \hat{\varphi}_i \hat{\varphi}_l \hat{\varphi}_j \,.
%> @f]
%>
%> @param  N    The local number of degrees of freedom, either as a scalar
%>              for all three basis functions, or as a vector with three
%>              entries, specifying the number of degrees of freedom for
%>              each basis function.
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D and gradPhi2D.
%> @param  qOrd (optional) The order of the quadrature rule to be used.
%> @retval ret  The computed array @f$[N\times N\times N\times 2]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = integrateRefElemDphiPhiPhi(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if length(N) == 1
  N = N * ones(3,1);
else
  validateattributes(N, {'numeric'}, {'numel', 3}, mfilename, 'N')
end % if
ret = zeros(N(1), N(2), N(3), 2);

if N(1) == 1 % p == 0
  return
end % if

if nargin < 3
  p = (sqrt(8*max(N)+1)-3)/2;  qOrd = max(2*p, 1); 
end % if
[~,~,W] = quadRule2D(qOrd);

for i = 1 : N(1)
  for j = 1 : N(2)
    for l = 1 : N(3)
      for m = 1 : 2
        ret(i,j,l,m) = sum( W' .* basesOnQuad.gradPhi2D{qOrd}(:,i,m) .* basesOnQuad.phi2D{qOrd}(:,j) .* basesOnQuad.phi2D{qOrd}(:,l) );
      end % for
    end % for
  end % for
end % for
end % function
