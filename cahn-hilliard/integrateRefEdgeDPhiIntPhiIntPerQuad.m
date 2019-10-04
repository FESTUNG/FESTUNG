% AUTOR: TIM ROITH
% Compute integrals over the edges of the reference element, whose integrands 
% consist of all permutations of two basis functions.

%===============================================================================
%> @file
%>
%> @brief Compute integrals over the edges of the reference element, whose 
%>        integrands consist of all permutations of two basis functions.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference element 
%>        @f$\hat{T}@f$, whose integrands consist of all permutations of two
%>        basis functions.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{{S}}}^\mathrm{diag}\in\mathbb{R}^{N\times N\times {n_\mathrm{edges}}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{S}}}^\mathrm{diag}]_{i,j,n} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_n(s) 
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
%> @retval ret  The computed array @f$[N\times N\times {n_\mathrm{edges}}]@f$
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
function ret = integrateRefEdgeDPhiIntPhiIntPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
nEdges = size(basesOnQuad.phi1D{end}, 3);

if nargin < 3
  switch nEdges
    case 3
      p = (sqrt(8*N+1)-3)/2;
    case 4
      p = sqrt(N)-1;
    otherwise
      error('Unknown number of edges')
  end % switch
  qOrd = 2*p+1;  
end % if

[~, W] = quadRule1D(qOrd); R = length(W);
ret = cell(2,1);
ret{1} = zeros(N, N, nEdges, R); % [N x N x nEdges x R]
ret{2} = zeros(N, N, nEdges, R); % [N x N x nEdges x R]

for n = 1 : nEdges
  for i = 1 : N
    for j = 1 : N
      for m = 1 : 2
        ret{m}(i,j,n,:) = ( basesOnQuad.gradPhi1D{qOrd}(:,i,n,m) .* basesOnQuad.phi1D{qOrd}(:,j,n) ) .* W.';
      end
    end  % for j
  end  % for i
end  % for n

end % function
