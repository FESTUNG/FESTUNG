% TODO
%===============================================================================
%> @file projectFuncCont2EdgeDataDisc.m
%>
%> @brief TODO
%===============================================================================
%>
%> @brief TODO
%>
%> TODO
%> 
%> All other entries are zero.
%> @param  g                TODO
%> 
%> @param  funcCont  TODO
%> 
%> @param  qOrd   TODO
%> 
%> @param  refEdgePhiPhi   TODO
%> 
%> @param  basesOnQuad   TODO
%> 
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
function dataDisc = projectFuncCont2EdgeDataDisc(g, funcCont, qOrd, refEdgePhiPhi, basesOnQuad)
validateattributes(funcCont, {'function_handle'}, {}, mfilename, 'funcCont');
validateattributes(basesOnQuad, {'struct'}, {}, mfilename, 'basesOnGamma');
qOrd = max(qOrd,1);  [Q, W] = quadRule1D(qOrd);
N = size(refEdgePhiPhi, 1);
rhs = zeros(g.numE, N);
for n = 1:3
  [Q1, Q2] = gammaMap(n, Q);
  rhs(g.E0T(:, n),:) = repmat(W, g.numT, 1) .* funcCont(g.mapRef2Phy(1, Q1, Q2), g.mapRef2Phy(2, Q1, Q2) ) * basesOnQuad.mu{qOrd};
end
dataDisc = rhs / refEdgePhiPhi;
end % function
