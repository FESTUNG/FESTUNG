% TODO
%===============================================================================
%> @file integrateRefElemDphiPhiPerQuad.m
%>
%> @brief NEW TODO
%===============================================================================
%>
%> @brief TODO
%>
%> TODO
%> 
%> All other entries are zero.
%> @param  N                TODO
%> 
%> @param  basesOnQuadEdge  TODO
%> 
%> @param  qOrd   TODO
%> 
%> @retval ret              TODO
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> 
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
%%TODO 
%
% I precompute matrices G_bar that allows for an 'easy' evaluation of 
% u_{m} phi_{kj} \partial_{x_m} phi_{ki}
%
function ret = integrateRefElemDphiPhiPerQuad(N, basesOnQuad, qOrd)
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnQuad')
if nargin < 3, p = (sqrt(8*N+1)-3)/2;  qOrd = max(2*p, 1);  end
[~, ~, W] = quadRule2D(qOrd); R = length(W);
ret = { zeros(N, N, R); zeros(N, N, R) };
for i = 1 : N
  for j = 1 : N
    for m = 1 : 2
      ret{m}(i, j, :) =  W(:) .* basesOnQuad.phi2D{qOrd}(:, j) .* basesOnQuad.gradPhi2D{qOrd}(:, i, m);
    end % for
  end % for
end % for
end % function
