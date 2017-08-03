% TODO

%===============================================================================
%> @file integrateRefEdgePhiIntMuPerQuad.m
%>
%> @brief NEW TODO
%===============================================================================
%>
%> @brief TODO
%>
%> 
%> All other entries are zero.
%> @param  g          The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param  markE0T    A marker indicating whether an edge should be 
%>                    recognized or not. @f$[K \times 3]@f$
%> @param  TODO TODO
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017
%> @author Balthasar Reuter, 2017
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
function ret = integrateRefEdgePhiIntMuPerQuad(N, basesOnQuad, qOrd)
validateattributes(N, {'numeric'}, {'numel', 2}, mfilename, 'N')
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = 2*p+1; end
[~, W] = quadRule1D(qOrd); R = length(W);
ret = { zeros(N(1), N(2), 3, R); zeros(N(1), N(2), 3, R) };
for n = 1 : 3 
  for i = 1 : N(1)
    for j = 1 : N(2)
      ret{1}(i, j, n, :) = W(:) .* basesOnQuad.phi1D{qOrd}(:, i, n) ...
                              .* basesOnQuad.mu{qOrd}(:, j);
      ret{2}(i, j, n, :) = W(:) .* basesOnQuad.phi1D{qOrd}(:, i, n) ...
                              .* basesOnQuad.thetaMu{qOrd}(:, j, 2);
    end % for j
  end % for i
end % for n
end % function
