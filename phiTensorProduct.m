% Evaluates the i-th basis function built from a tensor-product of two 
% one-dimensional basis functions.

%===============================================================================
%> @file phiTensorProduct.m
%>
%> @brief Evaluates the i-th basis function built from a tensor-product of two 
%>        one-dimensional basis functions.
%===============================================================================
%>
%> @brief Evaluates the i-th basis function built from a tensor-product of two 
%>        one-dimensional basis functions.
%>
%> @param  i          The index of the basis function.
%> @param  X1         A list of @f$\hat{x}^1@f$ coordinates.
%> @param  X2         A list of @f$\hat{x}^2@f$ coordinates.
%> @param  phiX1      Function handle to the one-dimensional basis function
%>                    in x1-direction.
%> @param  phiX2      Function handle to the one-dimensional basis function
%>                    in x2-direction.
%> @retval ret The @f$i@f$-th basis function in all points specified by
%>             <code>X1</code>, <code>X2</code>. It holds <code>size(X1) ==
%>             size(X2) == size(ret)</code>.
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
function [ ret ] = phiTensorProduct(i, X1, X2, phiX1, phiX2)
[m,n] = mapTensorProductIndexInv(i);
ret = phiX1(m, X1) .* phiX2(n, X2);
end  % function