% Evaluates the (spatial) derivatives of a basis function in all quadrature
% points on the reference triangle and multiplies with the according quadrature 
% weight.
%
%===============================================================================
%> @file integrateRefElemDphiPerQuad.m
%>
%> @brief Evaluates the (spatial) derivatives of a basis function in all 
%>				quadrature points on the reference triangle and multiplies with the 
%>				according quadrature weight.
%===============================================================================
%>
%> @brief Evaluates the (spatial) derivatives of a basis function in all 
%>				quadrature points on the reference triangle and multiplies with the 
%>				according quadrature weight.
%>
%> It computes a multidimensional array @f$\hat{\mathsf{F}}
%>    \in \mathbb{F}^{N\times R\times 2}@f$, which is defined by
%> @f[
%> \left[\hat{\mathsf{F}}\right]_{i,r,m} \;:=\;
%> \partial_{\hat{x}^m} \hat{\varphi}_i (\hat{q_r})\, (\hat{q_r})\, w_r
%> @f]
%> with the quadrature points @f$q_r@f$ given by
%> <code>quadRule2D()</code>
%>
%> @param  N    The local number of degrees of freedom
%> @retval ret  The computed array @f$[N\times R\times 2]@f$
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
function ret = integrateRefElemDphiPerQuad(N)
global gGradPhi2D
p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  [~, ~, W] = quadRule2D(qOrd);
ret = zeros(N, length(W), 2);
for i = 1 : N
	ret(i, :, 1) = gGradPhi2D{qOrd}(:, i, 1) .* W.';
	ret(i, :, 2) = gGradPhi2D{qOrd}(:, i, 2) .* W.';
end % for
end % function
