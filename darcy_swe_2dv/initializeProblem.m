% Calls the initialization routines of the subproblems.

%===============================================================================
%> @file
%>
%> @brief Calls the initialization routines of the subproblems.
%===============================================================================
%>
%> @brief Calls the initialization routines of the subproblems.
%>
%> This routine is called after darcy_swe_2dv/preprocessProblem.m.
%>
%> It executes darcy_2dv/initializeProblem.m and @link
%> swe_2dv/initializeProblem.m @endlink.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by configureProblem() and 
%>                      preprocessProblem(). @f$[\text{struct}]@f$
%>
%> @retval problemData  The input struct enriched with initial data.
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
function problemData = initializeProblem(problemData)
problemData.darcyData = problemData.darcySteps.initializeProblem(problemData.darcyData);
problemData.sweData = problemData.sweSteps.initializeProblem(problemData.sweData);
problemData.isFinished = false;
end % function