% Generates function handles for the steps in the problem solver framework.

%===============================================================================
%> @file getStepHandles.m
%>
%> @brief Generates function handles for the steps in the problem solver framework.
%===============================================================================
%>
%> @brief Generates function handles for the steps in the problem solver framework.
%>
%> It takes the specified `problemName` and generates function handles for
%> this problem for all the implemented steps.
%> The names of the steps are either provided as optional argument or given 
%> by getStepLists().
%>
%> @param  problemName      The name of the problem (i.e., name of the
%>                          folder with the implemented steps).
%> @param  stepList         (optional) The list of step names.
%>
%> @retval stepHandles      A struct containing the function handles for
%>                          each step
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
function stepHandles = getStepHandles(problemName, stepList)
if nargin < 2
  [preprocessList, stepList, postprocessList, subStepList] = getStepLists();
  if isequal(cellfun(@(fun) exist([problemName filesep fun '.m'], 'file'), subStepList), 2 * ones(size(subStepList)))
    stepList = [ preprocessList(:); stepList(:); subStepList(:); postprocessList(:) ];
  else
    stepList = [ preprocessList(:); stepList(:); postprocessList(:) ];
  end % if
end % if
stepHandles = struct();
for nFunc = 1 : length(stepList)
  stepHandles.(stepList{nFunc}) = getFunctionHandle([problemName filesep stepList{nFunc}]);
end % for
end % function