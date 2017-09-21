% Perturbes any kind of data with uniformally distributed noise of a certain 
% magnitude.

%===============================================================================
%> @file addNoise.m
%>
%> @brief Perturbes any kind of data with uniformally distributed noise of a 
%> 		  certain magnitude.
%===============================================================================
%>
%> @brief Perturbes any kind of data with uniformally distributed noise of a 
%> 		  certain magnitude.
%>
%> @param  data       Any kind of real-valued tensor, e.g. a representation of 
%> 					  a discrete function 
%>                    @f$d_h(\mathbf(x))@f$, e.g., as computed by 
%>                    <code>projectFuncCont2DataDisc()</code>
%>                    @f$[K \times N]@f$
%>
%> @param  absNoiseLvl The maximum possible value that perturbed data may differ
%>					   from original.
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
function data = addNoise(data, absNoiseLvl)
% data = (1+relNoiseLvl*2*(rand(size(data))-0.5)).*data;
data = absNoiseLvl*2*(rand(size(data))-0.5) + data;
end % function
