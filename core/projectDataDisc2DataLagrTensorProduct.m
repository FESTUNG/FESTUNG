% Convert the DG/modal into a Lagrange/nodal basis representation.

%===============================================================================
%> @file
%>
%> @brief Convert the DG/modal into a Lagrange/nodal basis representation.
%===============================================================================
%>
%> @brief Converts the DG/modal basis representation into a Lagrange/nodal 
%>        basis representation.
%>
%> @param  dataDisc  Coefficient matrix of the DG/modal basis representation.
%>                   @f$[K \times N]@f$
%> @param  p         (optional) Polynomial order in Lagrange basis (0, 1, or 2).
%>                   @f$[\text{scalar}]@f$
%> @retval The representation matrix of the Lagrange/nodal basis representation.
%>         @f$[K \times N]@f$ for @f$p \in \{0, 1, 2\}@f$ and
%>         @f$[K \times 9]@f$ for @f$p > 2@f$, i.e., all higher approximation
%>         orders are projected onto a quadratic Lagrange basis, since
%>         <code>visualizeDataLagrTrap()</code> can visualize up to elementwise
%>         quadratics only.
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
function dataLagr = projectDataDisc2DataLagrTensorProduct(dataDisc, p)
[K, N] = size(dataDisc);
if nargin == 1
  p = sqrt(N)-1;
end
switch p
  case 0,     L1 = 0.5;           L2 = 0.5;              % locally constant
  case 1,     L1 = [0, 1, 1, 0];  L2 = [0, 0, 1, 1];     % locally linear
  otherwise,  L1 = [0, 1, 1, 0, 0.5, 1, 0.5, 0, 0.5];    % locally quadratic
              L2 = [0, 0, 1, 1, 0, 0.5, 1, 0.5, 0.5];
end  % switch
dataLagr = zeros(K, length(L1));
for i = 1 : N
  dataLagr = dataLagr + dataDisc(:, i) * phiTensorProduct(i, L1, L2, @phi1D, @phi1D);
end  % for
end  % function