% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file darcyVert_sweVert/solveStep.m
%>
%> @brief Second step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed second in each loop iteration and is intended to
%> produce the solution at the next step, e.g., at a new time-level.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next step. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = solveStep(problemData, nStep)
%% Solve SWE time steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);

%% Solve Darcy time step
problemData.darcyData = problemData.darcySteps.preprocessStep(problemData.darcyData, nStep);

% Compute coupling terms
problemData.hSWE = reshape(1 / problemData.darcyData.tau * problemData.hSWE.', [], 1);
% for m = 1 : 2
%   problemData.darcyData.globJcouple{m} = problemData.darcyData.tildeGlobQcouple{m} * problemData.hSWE;
% end % for
% problemData.darcyData.globKcouple = problemData.darcyData.tildeGlobScouple * problemData.hSWE;

problemData.darcyData = problemData.darcySteps.solveStep(problemData.darcyData, nStep);
problemData.darcyData = problemData.darcySteps.postprocessStep(problemData.darcyData, nStep);
end % function