% Compute integrals over the edges of the reference element, whose integrands 
% consist of all permutations of three basis functions of which two are from
% a neighboring element.

%===============================================================================
%> @file
%>
%> @brief Compute integrals over the edges of the reference element, whose 
%>        integrands consist of all permutations of three basis functions of 
%>        which two are from a neighboring element.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference element 
%>        @f$\hat{T}@f$, whose integrands consist of all permutations of three
%>        basis functions, of which two belong to a neighboring element that
%>        are transformed using @f$\hat{\mathbf{\vartheta}}@f$.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{{R}}}^\mathrm{offdiag}
%> \in \mathbb{R}^{N\times N\times N\times {n_\mathrm{edges}}\times {n_\mathrm{edges}}}@f$, 
%> which is defined by
%> @f[
%> [\hat{\mathsf{{R}}}^\mathrm{offdiag}]_{i,j,l,n^-,n^+} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) 
%>   \hat{\varphi}_l\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s)
%>   \hat{\varphi}_j\circ \hat{\mathbf{\vartheta}}_{n^-n^+} \circ
%>   \hat{\mathbf{\gamma}}_{n^-}(s) \mathrm{d}s \,,
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code> and the mapping 
%> @f$\hat{\mathbf{\vartheta}}_{n^-n^+}@f$ as described in <code>theta()</code>.
%>
%> @param  N            The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D and thetaPhi1D.
%> @param  qOrd         (optional) The order of the quadrature rule to be used. 
%>                      Defaults to @f$2p+1 @f$.
%> @retval ret  The computed array 
%>              @f$[N\times N\times N \times {n_\mathrm{edges}}]@f$ or
%>              @f$[N\times N\times N\times {n_\mathrm{edges}}\times {n_\mathrm{edges}}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2017.
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
function ret = integrateRefEdgePhiIntPhiExtPhiExt(N, basesOnQuad, qOrd)
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

[~, W] = quadRule1D(qOrd); 

if nEdges == 3
  
  ret = zeros(N, N, N, 3, 3); % [N x N x N x 3 x 3]
  
  for nn = 1 : 3 % 3 edges
    for np = 1 : 3
      for l = 1 : N
        for i = 1 : N
          for j = 1 : N
            ret(i, j, l, nn, np) = sum( W'.* basesOnQuad.phi1D{qOrd}(:, i, nn) .* ...
              basesOnQuad.thetaPhi1D{qOrd}(:, l, nn, np) .* basesOnQuad.thetaPhi1D{qOrd}(:, j, nn, np) );
          end % for j
        end % for i
      end % for l
    end % for np
  end % for nn
  
elseif nEdges == 4
  
  ret = zeros(N, N, N, 4); % [N x N x N x 4]

  for nn = 1 : 4
    np = mapLocalEdgeIndexQuadri(nn);
    for l = 1 : N
      for i = 1 : N
        for j = 1 : N
          ret(i, j, l, nn) = sum( W' .* basesOnQuad.phi1D{qOrd}(:, i, nn) .* ...
            basesOnQuad.phi1D{qOrd}(:, l, np) .* basesOnQuad.phi1D{qOrd}(:, j, np) );
        end % for j
      end % for i
    end % for l
  end % for nn
  
end % if  
end % function
