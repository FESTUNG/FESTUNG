% Performs all post-processing steps, such as error estimates of the final
% solution and output of the final solution.

%===============================================================================
%> @file
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution and output of the final solution.
%===============================================================================
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution and output of the final solution.
%>
%> This routine is called after the main loop.
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
%% Visualization
if problemData.isVisSol
  cLagrange = projectDataDisc2DataLagr(problemData.cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, ...
                    problemData.numSteps, problemData.outputTypes);
end % if
%% Error evaluation
if isfield(problemData, 'cCont')
  problemData.error = computeL2Error(problemData.g, problemData.cDisc, ...
    @(x1, x2) problemData.cCont(problemData.tEnd, x1, x2), ...
    problemData.qOrd + 1, problemData.basesOnQuad);
  fprintf('L2 error w.r.t. the analytical solution: %g\n', problemData.error);
else
  problemData.error = computeL2Error(problemData.g, problemData.cDisc, ...
    problemData.c0Cont, problemData.qOrd + 1, problemData.basesOnQuad);
  fprintf('L2 error w.r.t. the initial condition: %g\n', problemData.error);
end % if
end % function

