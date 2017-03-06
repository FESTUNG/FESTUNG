% Evaluates a basis function in all quadrature points of the reference triangle
% and multiplies with the associated quadrature weight.

%===============================================================================
%> @file integrateRefElemPhiPerQuad.m
%>
%> @brief Evaluates a basis function in all quadrature points of the reference 
%>				triangle and multiplies with the associated quadrature weight.
%===============================================================================
%>
%> @brief Evaluates a basis function in all quadrature points of the reference 
%>				triangle and multiplies with the associated quadrature weight.
%>
%> It computes a matrix @f$\hat{\mathsf{P}}
%>    \in \mathbb{R}^{N \times R}@f$, which is defined by
%> @f[
%> \left[\hat{\mathsf{P}}\right]_{i,r} \;:=\;
%> \hat{\varphi}_i(\hat{q}_r) \hat{w}_r\,
%> @f]
%> with the quadrature points and weights @f$\hat{q}_r, \hat{w}_r@f$
%> of which there are @fR@f given by <code>quadRule2D()</code>
%>
%> @param  N            The local number of degrees of freedom
%> @param  basesOnQuad  A struct containing precomputed values of the basis
%>                      functions on quadrature points. Must provide at
%>                      least phi2D.
%> @retval ret  The computed array @f$[N\times R]@f$
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
function ret = integrateRefElemPhiPerQuad(N, basesOnQuad)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [~, ~, W] = execin('../quadRule2D',qOrd);
ret = zeros(N, length(W)); % [N x R]
for i = 1 : N
	ret(i, :) = basesOnQuad.phi2D{qOrd}(:,i) .* W.';
end % for
end % function
