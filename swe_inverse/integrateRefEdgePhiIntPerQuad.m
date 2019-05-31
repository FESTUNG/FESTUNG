% Evaluates a basis function in all quadrature points of all edges of the 
% reference triangle and multiplies with the according quadrature weight.

%===============================================================================
%> @file
%>
%> @brief Evaluates a basis function in all quadrature points of all edges of the 
%>				reference triangle and multiplies with the according quadrature weight.
%===============================================================================
%>
%> @brief Evaluates a basis function in all quadrature points of all edges of the 
%>				reference triangle and multiplies with the according quadrature weight.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{Q}}
%>    \in \mathbb{R}^{N\times R\times 3}@f$, which is
%> defined by
%> @f[
%> \left[\hat{\mathsf{Q}}\right]_{i,r,n} \;:=\;
%> \hat{\varphi}_i\circ\hat{\mathbf{\gamma}}_n(\hat{q}^r)\, \hat{w}^r
%> @f]
%> with the mapping @f$\hat{\mathbf{\gamma}}_n@f$ defined in 
%> <code>gammaMap()</code> and the quadrature points and weights @f$\hat{q}^r,
%> \hat{w}^r@f$ of which there are @fR@f given by <code>quadRule1D()</code>
%>
%> @param  N    The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi1D.
%> @retval ret  The computed array @f$[N\times R\times 3]@f$
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
function ret = integrateRefEdgePhiIntPerQuad(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p+1,1);  [~, W] = quadRule1D(qOrd);
ret = zeros(N, length(W), 3); % [N x R x 3]
for i = 1 : N
	for n = 1 : 3 % 3 edges
    ret(i,:,n) = basesOnQuad.phi1D{qOrd}(:,i,n) .* W.';
  end % for
end % for
end % function
