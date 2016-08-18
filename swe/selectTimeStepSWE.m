% Routine that adjusts the time step for the 2D Shallow-Water Equations.

%===============================================================================
%> @file selectTimeStepSWE.m
%>
%> @brief Routine that adjusts the time step for the 2D Shallow-Water Equations.
%===============================================================================
%>
%> @brief Routine that adjusts the time step for the 2D Shallow-Water Equations.
%>
%> This routine is specifically designed for the 2D Shallow-Water Equations to
%> compute a maximal time step size for which the model is stable. This is done
%> by the problems CFL condition. If the old time step is sufficient for 
%> stability it is increased to reduce the number of required time steps.
%>


%> TODO: input, output


%>
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, 
%>                      Vadym Aizinger
%>                      
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
function dt = selectTimeStepSWE(dx, dy, avgDepth, gConst, dt, nStep)
etal2 = 5; % TODO
 uul2 = 5; % TODO
 vvl2 = 5; % TODO
hav = avgDepth + etal2;
c = sqrt(gConst*hav);
umax = c + abs(uul2);
vmax = c + abs(vvl2);
dt1 = min(0.5 ./ (umax ./ dx + vmax ./ dy));
if dt1 < dt % TODO evtl andere Fallunterscheidung
	warning(['Time increment reduction necessary in step ' num2str(nStep) ': old value ' num2str(dt) ', new value = ' num2str(dt1) '.']);
end % if
dt = min([dt1 1.05*dt]);
end % function
