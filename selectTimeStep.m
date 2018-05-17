% Routine that adapts the time step for the 2D Shallow-Water Equations.

%===============================================================================
%> @file
%>
%> @brief Routine that adapts the time step for the 2D Shallow-Water Equations.
%===============================================================================
%>
%> @brief Routine that adapts the time step for the 2D Shallow-Water Equations.
%>
%> This routine computes the maximum time step size for which the model is 
%> stable using the CFL condition. If the old time step is sufficient for 
%> stability it is increased to reduce the number of required time steps.
%>
%> First the mean of the vertex values of the surface elevation and
%> velocities are computed for each triangle. Then the known mean depth is
%> used to compute @f$c = \sqrt{gH}@f$, where @f$g@f$ is the gravitaional
%> constant, and @f$H@f$ is the mean depth plus the mean surface elevation.
%> @f$c@f$ is added to the absolute value of velocity means, and these
%> values are divided by the mean edge dimensions for @f$x@f$-, and @f$y@f$
%> -components, respectively. Summed, inverted and weighted with @f$0.5@f$
%> to account for two-dimensions, one obtains an estimate on the maximum 
%> time step for each element. The minimum of these is the new time step.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct with the new time step size. 
%>                      @f$[\text{struct}]@f$
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
function problemData = selectTimeStep(problemData, nStep)
p = problemData.p;
qOrd = max(2*p,1);

hDisc = problemData.cDisc(:,:,1) - problemData.zbDisc;

meanElev = mean(projectDataDisc2DataLagr(problemData.cDisc(:,:,1), max(p,1)), 2);

dataQ0T = (problemData.cDisc(:,:,2) * problemData.basesOnQuad.phi2D{qOrd}.') ./ (hDisc * problemData.basesOnQuad.phi2D{qOrd}.');
dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
meanVelX = mean(projectDataDisc2DataLagr(dataDisc, max(p,1)), 2);

dataQ0T = (problemData.cDisc(:,:,2) * problemData.basesOnQuad.phi2D{qOrd}.') ./ (hDisc * problemData.basesOnQuad.phi2D{qOrd}.');
dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
meanVelY = mean(projectDataDisc2DataLagr(dataDisc, max(p,1)), 2);

hav = problemData.avgDepth + meanElev;
c = sqrt(problemData.gConst*hav);
umax = c + abs(meanVelX);
vmax = c + abs(meanVelY);
dt1 = min(0.5 ./ (umax ./ problemData.avgDiff(:,1) + vmax ./ problemData.avgDiff(:,2))) * 3 / (2*p+1);
if dt1 < problemData.dt
	warning(['Time increment reduction necessary in step ' num2str(nStep) ': old value ' num2str(problemData.dt) ', new value = ' num2str(dt1) '.']);
else
  warning(['Time increment increase performed in step ' num2str(nStep) ': old value ' num2str(problemData.dt) ', new value = ' num2str(dt1) '.']);
end % if
problemData.dt = dt1;
end % function
