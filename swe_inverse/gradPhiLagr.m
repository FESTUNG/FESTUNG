% Evaluates the m-th component of the i-th Lagrange basis function.
%
%===============================================================================
%> @file gradPhiLagr.m
%>
%> @brief Evaluates the m-th component of the i-th Lagrange basis function.
%===============================================================================
%>
%> @brief Evaluates the @f$m@f$-th component of the @f$i@f$-thLagrange basis 
%>        function on the reference triangle @f$\hat{T}@f$ at points specified
%>        by a list of @f$\hat{x}^1@f$ and @f$\hat{x}^2@f$ coordinates .
%>
%> @param  i   The index of the Lagrange basis function.
%> @param  m   The component of the gradient, @f$m\in\{1,2\}@f$.
%> @param  X1  A list of @f$\hat{x}^1@f$ coordinates.
%> @param  X2  A list of @f$\hat{x}^2@f$ coordinates.
%> @retval ret The @f$m@f$-th component of the @f$i@f$-th basis function in all 
%>             points specified by <code>X1</code>, <code>X2</code>. 
%>             It holds <code>size(X1) == size(X2) == size(ret)</code>.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Hennes Hajduk, Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = gradPhiLagr(i, m, X1, X2)
switch m
  case 1
    switch i
      case 1,  ret = -ones(size(X1));
      case 2,  ret = ones(size(X1));
      case 3,  ret = zeros(size(X1));
    end % switch
  case 2
    switch i
      case 1,  ret = -ones(size(X1));
      case 2,  ret = zeros(size(X1));
      case 3,  ret = ones(size(X1));
    end % switch
end % switch
end % function
