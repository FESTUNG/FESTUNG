% Second step of the four-part algorithm in the main loop. Carries out all
% Runge-Kutta steps for a time-step.

%===============================================================================
%> @file
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>        Carries out all Runge-Kutta steps for a time-step.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>        Carries out all Runge-Kutta steps for a time-step.
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
%> This routine obtains the parameters of the Runge-Kutta method and calls
%> iterateSubSteps() to actually carry out the Runge-Kutta method.
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
%> @author Alexander Jaust, 2017.
%> @author Balthasar Reuter, 2017.
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
if problemData.isStationary
  problemData.t = 0;
else
  % Obtain Runge-Kutta rule
  [problemData.t, problemData.A, problemData.b] = rungeKuttaImplicit(problemData.ordRK, problemData.dt, (nStep - 1) * problemData.dt);

  % Initialize solution vectors for RK steps
  problemData.cDiscRK = cell(length(problemData.t), 1);
end % if

% Carry out RK steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep, problemData.subStepHandles);
end % function