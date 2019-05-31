% First step of the three-part substepping algorithm.

%===============================================================================
%> @file
%>
%> @brief First step of the three-part substepping algorithm.
%===============================================================================
%>
%> @brief First step of the three-part substepping algorithm.
%>
%> The routine iterateSubSteps() repeatedly executes three steps until the 
%> parameter <code>problemData.isSubSteppingFinished</code> becomes 
%> <code>true</code>.
%> These three steps are:
%>
%>  1. darcy_swe_2dv/preprocessSubStep.m
%>  2. darcy_swe_2dv/solveSubStep.m
%>  3. darcy_swe_2dv/postprocessSubStep.m
%> 
%> This routine is executed first in each loop iteration and calls
%> @link swe_2dv/preprocessStep.m @endlink
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      @link darcy_swe_2dv/configureProblem.m @endlink and 
%>                      @link darcy_swe_2dv/preprocessProblem.m @endlink. 
%%>                     @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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
function problemData = preprocessSubStep(problemData, nStep, nSubStep)
problemData.sweData = problemData.sweSteps.preprocessStep(problemData.sweData, (nStep-1) * problemData.numSubSteps + nSubStep);
end % function
