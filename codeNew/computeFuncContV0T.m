% Assembles a matrix containing a function evaluated in vertex of each
% element.
%
%===============================================================================
%> @file computeFuncContV0T.m
%>
%> @brief Assembles a matrix containing a function evaluated in vertex of 
%>        each element.
%===============================================================================
%>
%> @brief Assembles a matrix containing a function evaluated in vertex of 
%>        each element.
%>
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  funcCont   A function handle for the continuous function.
%> @retval ret        The assembled matrix @f$[K \times 3]@f$
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
function valV0T = computeFuncContV0T(g, funcCont)
% Check function arguments that are directly used
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');

% Evaluate function
valV0T = zeros(g.numT,3);
for n = 1 : 3
    valV0T(:, n) = funcCont(g.coordV0T(:, n, 1), g.coordV0T(:, n, 2));
end
end