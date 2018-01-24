% Compute integrals over the edges of the reference element, whose integrands 
% consist of all permutations of three basis functions.

%===============================================================================
%> @file integrateRefEdgePhiIntPhiIntPhiInt.m
%>
%> @brief Compute integrals over the edges of the reference element, whose 
%>        integrands consist of all permutations of three basis functions.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference element 
%>        @f$\hat{T}@f$, whose integrands consist of all permutations of three
%>        basis functions.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{{R}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times N\times{n_\mathrm{edges}}}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{diag}]_{i,j,l,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_l \circ \hat{\mathbf{\gamma}}_n(s) 
%>   \hat{\varphi}_j\circ \hat{\mathbf{\gamma}}_n(s) \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMap()</code>.
%>
%> @param  N            The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @param  qOrd         (optional) The order of the quadrature rule to be used. 
%>                      Defaults to @f$2p+1 @f$.
%> @retval ret  The computed array @f$[N\times N\times N\times {n_\mathrm{edges}}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2017
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
function ret = integrateRefEdgePhiIntPhiIntPhiInt(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')

if length(N) == 1
  N = N * ones(3,1);
else
  validateattributes(N, {'numeric'}, {'numel', 3}, mfilename, 'N')
end % if

nEdges = size(basesOnQuad.phi1D{end}, 3);

if nargin < 3
  switch nEdges
    case 3
      p = (sqrt(8*max(N)+1)-3)/2;
    case 4
      p = sqrt(max(N))-1;
    otherwise
      error('Unknown number of edges')
  end % switch
  qOrd = 2*p+1;  
end % if

[~, W] = quadRule1D(qOrd);
ret = zeros(N(1), N(2), N(3), nEdges); % [N x N x N x nEdges]

for n = 1 : nEdges
  for l = 1 : N(3) % N basisfcts for D(t)
    for i = 1 : N(1)
      for j = 1 : N(2)
        ret(i,j,l,n) = sum(W' .* basesOnQuad.phi1D{qOrd}(:,i,n) .* basesOnQuad.phi1D{qOrd}(:,l,n) .* basesOnQuad.phi1D{qOrd}(:,j,n));
      end % for j
    end % for i
  end % for l
end % for n
end % function
