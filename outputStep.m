% Last step of the four-part algorithm in the main loop.

%===============================================================================
%> @file template/outputStep.m
%>
%> @brief Last step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Last step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the number of
%> iterations provided by configureProblem in the parameter
%> <code>numSteps</code> is reached. These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed last in each loop iteration and is intended to
%> provide output operations for the solution, e.g., write it to a file
%> for later visualization.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures, as provided 
%>                      by configureProblem() and preprocessProblem(). 
%>                      @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with post-processed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = outputStep(problemData, nStep)
%% Visualize solution and evaluate stations.
problemData = problemData.visualizeSolution(problemData, nStep);

%% Hot-start file output
if problemData.isHotstartOutput && mod(nStep, problemData.hotstartOutputFrequency) == 0
  hotstartStep = mod(nStep / problemData.hotstartOutputFrequency, 2);
  outputHotstart(['output/' problemData.name '_' num2str(hotstartStep)], 'cDisc', problemData.cDisc, 't', problemData.t);
end % if

%% Update waitbar.
if problemData.isWaitbar
  percentDone = round( nStep / problemData.numSteps * 100 );
  problemData.waitbar = waitbar(percentDone / 100, problemData.waitbar, strcat( [ 'Time stepping:', ' ', num2str(percent), str ] ) );
end % if
end % function

