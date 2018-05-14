% Assembles a multidimensional array containing a function evaluated in
% each quadrature point of each edge, multiplied with the edge normal.

%===============================================================================
%> @file
%>
%> @brief Assembles a multidimensional array containing a function evaluated
%>        in each quadrature point of each edge, multiplied with the edge normal.
%===============================================================================
%>
%> @brief Assembles a multidimensional array containing a function evaluated
%>        in each quadrature point of each edge, multiplied with the edge normal.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  funcCont1  A function handle for the @f$x^1@f$-component of the
%>                    continuous function.
%> @param  funcCont2  A function handle for the @f$x^2@f$-component of the
%>                    continuous function.
%> @param  qOrd       The order of the 1D-quadrature rule to be used. 
%>                    Determines number and position of quadrature points.
%> @retval ret        The assembled array @f$[K \times 3 \times R]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function ret = computeFuncContNuOnQuadEdge(g, funcCont1, funcCont2, qOrd)
% Extract dimensions and determine quadrature rule
K = g.numT;  [Q, W] = quadRule1D(qOrd);

% Check function arguments that are directly used
validateattributes(funcCont1, {'function_handle'}, {}, mfilename, 'funcCont1');
validateattributes(funcCont2, {'function_handle'}, {}, mfilename, 'funcCont2');

% Evaluate function
ret = zeros(K, 3, length(W));
for n = 1 : 3
  [Q1, Q2] = gammaMap(n, Q);
  ret(:, n, :) = bsxfun(@times, g.nuE0T(:, n, 1), funcCont1(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2))) + ...
                 bsxfun(@times, g.nuE0T(:, n, 2), funcCont2(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2)));
end % for
end