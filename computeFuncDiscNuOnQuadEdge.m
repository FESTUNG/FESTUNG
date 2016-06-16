% Assembles a multidimensional array containing the values of a discrete function 
% in each quadrature point of each edge, multiplied with the edge normal.
%
%===============================================================================
%> @file computeFuncDiscNuOnQuadEdge.m
%>
%> @brief Assembles a multidimensional array containing the values of a discrete 
%>        function in each quadrature point of each edge, multiplied with the edge normal.
%===============================================================================
%>
%> @brief Assembles a multidimensional array containing the values of a discrete 
%>        function in each quadrature point of each edge, multiplied with the edge normal.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  dataDisc1  A representation of the @f$x^1@f$-component of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N_\mathrm{data}]@f$ 
%> @param  dataDisc2  A representation of the @f$x^2@f$-component of the discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N_\mathrm{data}]@f$ 
%> @param  qOrd       The order of the 1D-quadrature rule to be used. 
%>                    Determines number and position of quadrature points.
%> @retval ret        The assembled array @f$[K \times 3 \times R]@f$
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
function ret = computeFuncDiscNuOnQuadEdge(g, dataDisc1, dataDisc2, qOrd)
% Extract dimensions and determine quadrature rule
K = g.numT;  N = size(dataDisc1, 2);  [Q, W] = quadRule1D(qOrd);

% Check function arguments that are directly used
validateattributes(dataDisc1, {'numeric'}, {'size', [K N]}, mfilename, 'dataDisc1');
validateattributes(dataDisc2, {'numeric'}, {'size', [K N]}, mfilename, 'dataDisc2');

% Evaluate function
ret = zeros(K, 3, length(W));

for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  dataLagr = zeros(2*K, length(W));
  for i = 1 : N
    dataLagr = dataLagr + [dataDisc1(:, i); dataDisc2(:, i)] * phi(i, Q1, Q2);
  end % for
  ret(:, n, :) = bsxfun(@times, g.nuE0T(:, n, 1), dataLagr(1:K,:)) + ...
                 bsxfun(@times, g.nuE0T(:, n, 2), dataLagr(K+1:2*K,:));
end % for
end % function