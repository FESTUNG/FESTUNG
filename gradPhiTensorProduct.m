% Evaluates the m-th component of the gradient of the i-th basis function.

%===============================================================================
%> @file gradPhiTensorProduct.m
%>
%> @brief Evaluates the m-th component of the gradient of the i-th basis 
%>        function.
%===============================================================================
%>
%> @brief Evaluates the @f$m@f$-th component of the gradient of the @f$i@f$-th
%>        basis function, which is built by a tensor product from two one-
%>        dimensional basis functions, in points located in the respective 
%>        reference intervals of the one-dimensional basis functions, 
%>        specified by a list of @f$\hat{x}^1@f$ and @f$\hat{x}^2@f$ coordinates.
%>
%> @param  i          The index of the basis function.
%> @param  m          The component of the gradient, @f$m\in\{1,2\}@f$.
%> @param  X1         A list of @f$\hat{x}^1@f$ coordinates.
%> @param  X2         A list of @f$\hat{x}^2@f$ coordinates.
%> @param  phiX1      Function handle to the one-dimensional basis function
%>                    in x1-direction.
%> @param  phiX2      Function handle to the one-dimensional basis function
%>                    in x2-direction.
%> @param  gradPhiX1  Function handle to the gradient of the
%>                    one-dimensional basis function in x1-direction.
%> @param  gradPhiX2  Function handle to the gradient of the
%>                    one-dimensional basis function in x2-direction.
%> @retval ret The @f$m@f$-th component of the gradient of the @f$i@f$-th 
%>             basis function in all points specified by <code>X1</code>, 
%>             <code>X2</code>. 
%>             It holds <code>size(X1) == size(X2) == size(ret)</code>.
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
function ret = gradPhiTensorProduct(i, m, X1, X2, phiX1, phiX2, gradPhiX1, gradPhiX2)
[iX, iY] = mapTensorProductIndexInv(i);
switch m
  case 1
    ret = gradPhiX1(iX, X1) .* phiX2(iY, X2);
  case 2
    ret = phiX1(iX, X1) .* gradPhiX2(iY, X2);
end % switch m
end  % function

