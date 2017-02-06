% Fills the problem's data structures with initial data.

%===============================================================================
%> @file advection/initializeProblem.m
%>
%> @brief Fills the problem's data structures with initial data.
%===============================================================================
%>
%> @brief Fills the problem's data structures with initial data.
%>
%> This routine is called after advection/preprocessProblem.m.
%>
%> Before entering the main loop of the solution algorithm, this routine
%> fills the problem's data structures with initial data.
%>
%> It projects the initial data and performs slope limiting on the
%> projected function.
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
function problemData = initializeProblem(problemData)
problemData.isFinished = false;
%% Initial data on elements
problemData.cDisc = projectFuncCont2DataDisc(problemData.g, problemData.c0Cont, 2*problemData.p+1, ...
                                             problemData.hatM, problemData.basesOnQuad);
fprintf('L2 error w.r.t. the initial condition: %g\n', ...
  computeL2Error(problemData.g, problemData.cDisc, problemData.c0Cont, 2*problemData.p, problemData.basesOnQuad));

%% Initial data on faces
problemData.lamDisc = projectFuncCont2FaceDataDisc(problemData.g, problemData.c0Cont, 2*problemData.p+1, ...
                                                   problemData.hatMlambda, problemData.basesOnGamma);
% fprintf('L2 error w.r.t. the initial condition: %g\n', ...
%   computeL2Error(problemData.g, problemData.lamDisc, problemData.c0Cont, 2*problemData.p, problemData.basesOnQuad));

%% visualization of inital condition.
if problemData.isVisSol
  cLagrange = projectDataDisc2DataLagr(problemData.cDisc);
  visualizeDataLagr(problemData.g, cLagrange, 'u_h', problemData.outputBasename, 0, problemData.outputTypes)
end
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', ...
  problemData.tEnd, problemData.tau, problemData.numSteps)
end % function