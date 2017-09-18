% Mapping from interval [0,1] to the n-th edge in the reference triangle.

%===============================================================================
%> @file gammaMap.m
%>
%> @brief Mapping from interval [0,1] to the n-th edge in the reference triangle.
%===============================================================================
%>
%> @brief Computes the reference coordinates
%>        @f$\hat{\mathbf{x}} = (\hat{x}^1,\hat{x}^2)^T@f$ as given by the 
%>        mapping @f$\hat{\mathbf{\gamma}}_n(s)@f$.
%>
%> The mapping @f$\hat{\mathbf{\gamma}}_n : [0,1] \mapsto \hat{E}_n@f$ 
%> maps the unit interval to the @f$n@f$-th edge of the reference triangle
%> @f$\hat{T} = \{ (0,0), (1,0), (0,1) \}@f$ and is defined as
%> @f[
%> \hat{\mathbf{\gamma}}_1(s) = (1-s, s)^T \,,\quad
%> \hat{\mathbf{\gamma}}_2(s) = (0, 1-s)^T \,,\quad
%> \hat{\mathbf{\gamma}}_3(s) = (s,   0)^T \,.
%> @f]
%>
%> @param  n  The index of the edge in the reference triangle @f$\hat{T}@f$.
%>            @f$[\text{scalar}]@f$
%> @param  S  The parameter @f$s@f$ of the mapping. Can be a vector, e.g.,
%>            to compute the mapping for multiple values in one function call.
%> @retval X1,X2  reference coordinates
%>                @f$\hat{\mathbf{x}} = (\hat{x}^1,\hat{x}^2)^T@f$.
%>                It holds <code>size(X1) == size(X2) == size(S)</code>.
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
function [X1, X2] = gammaMap(n, S)
S = S(:)';
switch n
  case 1,  X1 = 1 - S;           X2 = S;
  case 2,  X1 = zeros(size(S));  X2 = 1 - S;
  case 3,  X1 = S;               X2 = zeros(size(S));
end % switch
end % function
