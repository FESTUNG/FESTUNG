% Performs all post-processing steps, such as error estimates of the final
% solution, etc.

%===============================================================================
%> @file
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%===============================================================================
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%>
%> This routine is called after the main loop.
%>
%> It can include things such as error estimates, an output operation of
%> the final solution, etc.
%>
%> @param  problemData  A struct with problem parameters and solution
%>                      vectors. @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description and precomputed fields.
%>                      @f$[\text{struct}]@f$
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
function problemData = postprocessProblem(problemData)
problemData.darcyData = problemData.darcySteps.postprocessProblem(problemData.darcyData);

if problemData.isCouplingSWE
  problemData = execin([problemData.problemName filesep 'preprocessStep'], problemData, problemData.numSteps + 1);
end % if
 
problemData.sweData = problemData.sweSteps.postprocessProblem(problemData.sweData);
if isfield(problemData.sweData, 'error') && isfield(problemData.darcyData, 'error')
  problemData.error = [ problemData.sweData.error, problemData.darcyData.error ];
elseif isfield(problemData.sweData, 'error')
  problemData.error = [ problemData.sweData.error, 0, 0, 0 ];
elseif isfield(problemData.darcyData, 'error')
  problemData.error = [ 0, 0, 0, problemData.darcyData.error ];
end % if
end

