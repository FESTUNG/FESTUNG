% Third step of the three-part substepping algorithm.

%===============================================================================
%> @file darcy_swe_2dv/postprocessSubStep.m
%>
%> @brief Third step of the three-part substepping algorithm.
%===============================================================================
%>
%> @brief Third step of the three-part substepping algorithm.
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
%> This routine is executed third in each loop iteration and is intended to
%> post-process the solution computed by solveSubStep().
%> It is also commonly the right place to decide whether the substepping
%> method has finished and to update
%> <code>problemData.isSubSteppingFinished</code> accordingly.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with postprocessed data
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
function problemData = postprocessSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.sweSteps.postprocessStep(problemData.sweData, (nStep - 1) * problemData.numSubSteps + nSubStep);
problemData.sweData = problemData.sweSteps.outputStep(problemData.sweData, (nStep - 1) * problemData.numSubSteps + nSubStep);

if problemData.isCouplingDarcy
  % Integrate water height over time for coupling (using trapezoidal rule)
  problemData.hSWE = problemData.hSWE + 0.5 * problemData.sweData.tau * (problemData.sweData.cDiscRK{1, 1} + problemData.sweData.cDiscRK{end, 1});
end % if

problemData.isSubSteppingFinished = nSubStep >= problemData.numSubSteps;
end % function
