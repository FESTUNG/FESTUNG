% Third step of the three-part substepping algorithm for each Runge-Kutta stage.

%===============================================================================
%> @file
%>
%> @brief Third step of the three-part substepping algorithm for each Runge-Kutta stage.
%===============================================================================
%>
%> @brief Third step of the three-part substepping algorithm for each Runge-Kutta stage.
%>
%> The routine iterateSubSteps() repeatedly executes three steps until the 
%> parameter <tt>problemData.isSubSteppingFinished</tt> becomes 
%> <tt>true</tt>.
%> These three steps are:
%>
%>  1. @link swe_2dv/preprocessSubStep.m @endlink
%>  2. @link swe_2dv/solveSubStep.m @endlink
%>  3. @link swe_2dv/postprocessSubStep.m @endlink
%> 
%> This routine is executed third in each loop iteration.
%> It decides whether the substepping is finished and updates
%> <tt>problemData.isSubSteppingFinished</tt> accordingly.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      @link swe_2dv/configureProblem.m @endlink and 
%>                      @link swe_2dv/preprocessProblem.m @endlink. 
%%>                     @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with postprocessed data
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
function problemData = postprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
problemData.isSubSteppingFinished = nSubStep >= length(problemData.omega);
end % function
