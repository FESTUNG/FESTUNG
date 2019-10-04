% Evaluates basis functions in the corners of the reference triangle and
% stores them in a struct.

%===============================================================================
%> @file computeBasesOnQuad.m
%>
%> @brief Evaluates basis functions in the corners of the reference triangle and
%>        stores them in a struct.
%===============================================================================
%>
%> @brief Evaluates basis functions in the corners of the reference triangle and
%>        stores them in a struct.
%>
%> It evaluates the basis functions provided by <code>phi()</code> in all corners 
%> on the reference triangle @f$\hat{T} = \{(0,0), (1,0), (0,1) \}@f$
%> 
%> This function computes the following struct variables (dimensions given for
%> each order):
%> - <code>phi2D</code>: @f$\hat{\varphi}_i(\mathbf{q}_r) \; [R \times N]@f$
%> 
%> @param  N          The number of local degrees of freedom. For polynomial
%>                    order @f$p@f$, it is given as @f$N = (p+1)(p+2)/2@f$
%>                    @f$[\text{scalar}]@f$
%>
%> @retval  basesOnQuad A struct with the computed array.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2015 Florian Frank, Balthasar Reuter, Vadym Aizinger
%>
%> @author Andreas Rupp, 2018.
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
function phiOnVertices = computeBasesOnCorners(N)

phiOnVertices = zeros(N,3);

for i = 1: N
  phiOnVertices(i,1) = phi(i,0,0);
  phiOnVertices(i,2) = phi(i,1,0);
  phiOnVertices(i,3) = phi(i,0,1);
end

end % function computeBasesOnCorners
