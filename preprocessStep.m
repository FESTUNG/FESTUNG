% First step of the four-part algorithm in the time stepping loop.

%===============================================================================
%> @file
%>
%> @brief First step of the four-part algorithm in the time stepping loop. 
%===============================================================================
%>
%> @brief First step of the four-part algorithm in the time stepping loop. 
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. @link swe_2dv/preprocessStep.m @endlink
%>  2. @link swe_2dv/solveStep.m @endlink
%>  3. @link swe_2dv/postprocessStep.m @endlink
%>  4. @link swe_2dv/outputStep.m @endlink
%> 
%> This routine is executed first in each loop iteration.
%> It takes triggers the adaptation of the free surface (see swe_2dv/adaptFreeSurface.m)
%> and initializes the solution vectors for the Runge-Kutta method.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      @link swe_2dv/configureProblem.m @endlink and 
%>                      @link swe_2dv/preprocessProblem.m @endlink. 
%%>                     @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct without any modifications.
%>                      @f$[\text{struct}]@f$
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
function problemData = preprocessStep(problemData, nStep)
% Apply mesh adaptation to free surface movement
problemData = problemData.fn_adaptFreeSurface(problemData);

% Copy last time step
problemData.cDiscRK(1, :) = problemData.cDiscRK(end, :);
end % function