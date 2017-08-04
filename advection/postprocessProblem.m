% Performs all post-processing steps, such as error estimates of the final
% solution and output of the final solution.

%===============================================================================
%> @file ./advection/postprocessProblem.m
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
function problemData = postprocessProblem(problemData)
%% Visualization
if problemData.isVisSol
  cLagrange = projectDataDisc2DataLagr(problemData.cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, ...
                    problemData.numSteps, problemData.outputTypes);
end % if
%% Error evaluation
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2Error(problemData.g, problemData.cDisc, problemData.c0Cont, 2*problemData.p, problemData.basesOnQuad));
fprintf('norm(cDisc, 1) = %g\n', norm(problemData.cDisc(:), 1));
fprintf('norm(cDisc, 2) = %g\n', norm(problemData.cDisc(:), 2));
fprintf('norm(cDisc, inf) = %g\n', norm(problemData.cDisc(:), inf));


end % function

