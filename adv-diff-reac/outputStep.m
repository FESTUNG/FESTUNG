% Last step of the four-part algorithm in the main loop.

%===============================================================================
%> @file
%>
%> @brief Last step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Last step of the four-part algorithm in the main loop.
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
%> This routine is executed last in each loop iteration and is intended to
%> provide output operations for the solution, e.g., write it to a file
%> for later visualization.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures, as provided 
%>                      by configureProblem() and preprocessProblem(). 
%>                      @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with post-processed data
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
function problemData = outputStep(problemData, nStep)
if problemData.isVisSol && mod(nStep, problemData.outputFrequency) == 0
  uLagr = projectDataDisc2DataLagr(problemData.uDisc);
  visualizeDataLagr(problemData.g, uLagr, 'u_h', problemData.outputBasename, ...
                    nStep / problemData.outputFrequency, problemData.outputTypes);
end % if
end % function

