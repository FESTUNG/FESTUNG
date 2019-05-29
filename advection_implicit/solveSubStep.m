% Compute the solution of the current Runge-Kutta stage.

%===============================================================================
%> @file
%>
%> @brief Compute the solution of the current Runge-Kutta stage.
%===============================================================================
%>
%> @brief Compute the solution of the current Runge-Kutta stage.
%>
%> The routine iterateSubSteps() repeatedly executes three steps until the 
%> parameter <code>problemData.isSubSteppingFinished</code> becomes 
%> <code>true</code>.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%> 
%> This routine is executed second in each loop iteration.
%> It assembles the global system and computes the solution at
%> the next Runge-Kutta level.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with the new solution
%>                      for this Runge-Kutta stage. @f$[\text{struct}]@f$
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
function problemData = solveSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
% Building the system
sysA = -problemData.globG{1} - problemData.globG{2} + problemData.globR;
sysV = problemData.globL - problemData.globKD - problemData.globKN;

if problemData.isStationary
  problemData.cDiscRK{nSubStep + 1} = sysA \ sysV;
else
  % Computing the rhs
  sysR = (problemData.tau * problemData.A(nSubStep, nSubStep)) * sysV + problemData.globM * problemData.cDiscRK{1};
  for j = 1 : nSubStep - 1
    sysR = sysR + (problemData.tau * problemData.A(nSubStep, j)) * problemData.rhsRK{j};
  end % for

  % Compute the next step
  problemData.cDiscRK{nSubStep + 1} = (problemData.globM + (problemData.tau * problemData.A(nSubStep, nSubStep)) * sysA) \ sysR; 

  % Store the new rhs
  problemData.rhsRK{nSubStep} = sysV - sysA * problemData.cDiscRK{nSubStep + 1};
end % if
end % function
