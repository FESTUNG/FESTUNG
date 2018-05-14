% Compute integrals over the edges of the reference square, whose integrands 
% consist of all permutations of a two-dimensional basis function and a 
% two-dimensional and a one-dimensional basis function from the
% neighbouring element.

%===============================================================================
%> @file
%>
%> @brief Compute integrals over the edges of the reference square, whose 
%>        integrands consist of all permutations of a two-dimensional basis 
%>        function and a two-dimensional and a one-dimensional basis function 
%>        from the neighbouring element.
%===============================================================================
%>
%> @brief Compute integrals over the edges of the reference square, whose 
%>        integrands consist of all permutations of a two-dimensional basis 
%>        function and a two-dimensional and a one-dimensional basis function 
%>        from the neighbouring element.
%>
%> It computes a multidimensional array
%> @f$\hat{\mathsf{{P}}}^\mathrm{offdiag}\in\mathbb{R}^{N\times\overline{N}\times4}@f$
%> defined by
%> @f[
%> [\hat{\mathsf{{P}}}^\mathrm{offdiag}]_{i,j,l,n^-} =
%>   \int_0^1 \hat{\varphi}_i \circ \hat{\mathbf{\gamma}}_{n^-}(s) \,
%>   \hat{\varphi}_j \circ \hat{\mathbf{\gamma}}_{n^+}(s) 
%>   \hat{\phi}_l\circ [\hat{\mathbf{\gamma}}_{n^+}(s)]_1 \mathrm{d}s \,,
%> @f]
%> where the mapping @f$\hat{\mathbf{\gamma}}_n@f$ is given in 
%> <code>gammaMapTetra()</code> and the index @f$n^+@f$ is given implicitely
%> as described in <code>mapLocalEdgeTetra()</code>.
%>
%> @param  N            The local number of degrees of freedom 
%>                      @f$\mathbf{N} = [N, \overline{N}]@f$.
%> @param  qOrd         The order of the quadrature rule to be used.
%> @param  basesOnQuad2D A struct containing precomputed values of the 
%>                      two-dimensional basis functions on quadrature points.
%>                       Must provide at least phi1D.
%> @param  basesOnQuad1D A struct containing precomputed values of the 
%>                       one-dimensional basis functions on quadrature points.
%>                       Must provide at least phi1D and phi0D.
%> @retval ret  The computed array @f$[N\times \overline{N}\times 4]@f$
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
function ret = integrateRefEdgeTetraPhiIntPhiExtPhi1DExt(N, qOrd, basesOnQuad2D, basesOnQuad1D)
[~, W] = quadRule1D(qOrd);
mapE0E = [2 1 4 3];
ret = zeros(N(1), N(1), N(2), 4);
for n = 1 : 2
  np = mapE0E(n);
  for i = 1 : N(1)
    for j = 1 : N(1)
      for l = 1 : N(2)
        ret(i,j,l,n) = W * ( basesOnQuad2D.phi1D{qOrd}(:,i,n) .* basesOnQuad2D.phi1D{qOrd}(:,j,np) .* basesOnQuad1D.phi1D{qOrd}(:,l) );
      end % for l
    end  % for j
  end  % for i
end  % for n
for n = 3 : 4
  np = mapE0E(n);
  for i = 1 : N(1)
    for j = 1 : N(1)
      for l = 1 : N(2)
        ret(i,j,l,n) = W * ( basesOnQuad2D.phi1D{qOrd}(:,i,n) .* basesOnQuad2D.phi1D{qOrd}(:,j,np) * basesOnQuad1D.phi0D{qOrd}(l,n-2) );
      end % for l
    end  % for j
  end  % for i
end  % for n
end  % function