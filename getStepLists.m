% Obtain the list of step names in the generic problem framework.

%===============================================================================
%> @file getStepLists.m
%>
%> @brief Obtain the list of step names in the generic problem framework.
%===============================================================================
%>
%> @brief Obtain the list of step names in the generic problem framework.
%>
%> It provides the names for all steps in the 
%> [generic problem framework](doxygen/solver-framework).
%>
%> @retval preprocessList   A cell array containing the steps in the
%>                          preprocessing phase.
%> @retval stepList         A cell array containing the steps in the
%>                          iterative solver phase.
%> @retval postprobessList  A cell array containing the steps in the
%>                          postprocessing phase.
%> @retval subStepList      A cell array containing the steps in the
%>                          optional substepping phase.
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Florian Frank, Balthasar Reuter, Vadym Aizinger
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
function [preprocessList, stepList, postprocessList, subStepList] = getStepLists()
preprocessList = { 'configureProblem'; 'preprocessProblem'; 'initializeProblem' };
stepList = { 'preprocessStep'; 'solveStep'; 'postprocessStep'; 'outputStep' };
postprocessList = { 'postprocessProblem' };
subStepList = { 'preprocessSubStep'; 'solveSubStep'; 'postprocessSubStep' };
end % function