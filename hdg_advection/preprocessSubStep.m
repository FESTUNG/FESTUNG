% Preprocessing of the Runge-Kutta step.

%===============================================================================
%> @file advection/preprocessSubStep.m
%>
%> @brief Preprocessing of the Runge-Kutta step.
%===============================================================================
%>
%> @brief Preprocessing of the Runge-Kutta step.
%>
%> The routine iterateSubSteps() repeatedly executes three steps until the 
%> parameter <code>problemData.isSubSteppingFinished</code> becomes 
%> <code>true</code>.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%> 
%> This routine is executed first in each loop iteration.
%> It takes care of the assembly of time-dependent matrices and right hand
%> side vectors. Furthermore, it evaluates the boundary conditions for the 
%> current Runge-Kutta step.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with preprocessed data
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
function problemData = preprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
K = problemData.K;
N = problemData.N;

% L2 projections of algebraic coefficients
fDisc  = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(problemData.t(nSubStep),x1,x2),  ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(problemData.t(nSubStep),x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(problemData.t(nSubStep),x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);

% Evaluate normal velocity in quadrature points of edges
vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(problemData.t(nSubStep),x1,x2), ...
                      @(x1,x2) problemData.u2Cont(problemData.t(nSubStep),x1,x2), 2*problemData.p+1);

% Assembly of time-dependent global matrices
problemData.globG = assembleMatElemDphiPhiFuncDiscVec(problemData.g, problemData.hatG, u1Disc, u2Disc);
problemData.globR = assembleMatEdgePhiPhiValUpwind(problemData.g, ~problemData.g.markE0TbdrN, ...
                                                   problemData.hatRdiagOnQuad, problemData.hatRoffdiagOnQuad, ...
                                                   vNormalOnQuadEdge, problemData.g.areaE0TbdrNotN);
                                               
% Assembly of Dirichlet boundary contributions
problemData.globKD = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
                        @(x1,x2) problemData.cDCont(problemData.t(nSubStep),x1,x2), vNormalOnQuadEdge, N, ...
                        problemData.basesOnQuad, problemData.g.areaE0TbdrD);

% Assembly of Neumann boundary contributions
gNUpwind = @(x1,x2) (problemData.gNCont(problemData.t(nSubStep),x1,x2) <= 0) .* problemData.gNCont(problemData.t(nSubStep),x1,x2);
problemData.globKN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                                                   gNUpwind, N, problemData.basesOnQuad);

% Assembly of the source contribution
problemData.globL = problemData.globM * reshape(fDisc', K*N, 1);

%% HDG stuff



end % function
