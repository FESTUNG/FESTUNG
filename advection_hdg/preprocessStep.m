% First step of the four-part algorithm in the main loop. Applies the mass
% matrix to the solution at the old time level.

%===============================================================================
%> @file
%>
%> @brief First step of the four-part algorithm in the main loop. 
%>        Applies the mass matrix to the solution at the old time level.
%===============================================================================
%>
%> @brief First step of the four-part algorithm in the main loop.
%>        Applies the mass matrix to the solution at the old time level.
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
%> This routine is executed first in each loop iteration.
%> The Advection problem requires substepping due to the Runge-Kutta method
%> (see solveStep() and @ref JRASK2018 for details). Thus, no terms can be
%> assembled here. However, every Runge-Kutta stage requires the solution
%> at the old time level with the mass matrix applied to it, which is
%> computed here.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct without any modifications.
%>                      @f$[\text{struct}]@f$
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
function problemData = preprocessStep(problemData, nStep) %#ok<INUSD>
if ~problemData.isStationary
  problemData.globMcDisc = problemData.globMphi * reshape(problemData.cDisc', [], 1);
end % if
end % function
