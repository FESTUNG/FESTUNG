% Perturbes a DG function with uniformally distributed noise of a certain 
% magnitude.
%
%===============================================================================
%> @file addNoiseInVertices2DataDisc.m
%>
%> @brief Perturbes a DG function with uniformally distributed noise of a 
%>        certain magnitude.
%===============================================================================
%>
%> @brief Perturbes a DG function with uniformally distributed noise of a 
%>        certain magnitude.
%>
%> @param g           The lists describing the geometric and topological 
%>                    properties of a triangulation (see 
%>                    <code>generateGridData()</code>) 
%>                    @f$[1 \times 1 \text{ struct}]@f$
%> @param dataDisc    A representation of the discrete function ,e.g., as 
%>                    computed by <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times {N}]@f$
%> @param absNoiseLvl The maximum possible value that perturbed data may 
%>                    differ from original.
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
function data = addNoiseInVertices2DataDisc(g, data, absNoiseLvl)
noiseV = absNoiseLvl*2*(rand(g.numV, 1)-0.5);
data = data + noiseV(g.V0T);
end % function