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
function problemData = preprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSD>



K = problemData.K;
N = problemData.N;


% Evaluate normal advection velocity on every edge.
uNormalQ0E0T = computeFuncContNuOnQuadEdge(problemData.g, problemData.u1Cont, problemData.u2Cont, ...
                  2*problemData.p+1);

% Assemble source term.
problemData.globH = problemData.globMphi * reshape(projectFuncCont2DataDisc(problemData.g, problemData.fCont, ...
                        max(2*problemData.p, 1), problemData.hatM, problemData.basesOnQuad)', [], 1);

% Assemble element integral contributions.
problemData.globG = assembleMatElemDphiPhiFuncContVec(problemData.g, problemData.hatG, ...
                        problemData.u1Cont, problemData.u2Cont);

% Assemble flux on interior edges
problemData.globS = assembleMatEdgePhiIntMuVal(problemData.g, problemData.g.markE0Tint, ...
                        problemData.hatS, uNormalQ0E0T);

% Assemble Dirichlet boundary conditions.
cQ0E0T = computeFuncContOnQuadEdge(problemData.g, problemData.cDCont, 2 * problemData.p + 1);
problemData.globFphiD = assembleVecEdgePhiIntVal(problemData.g, problemData.g.markE0TbdrD, ...
                            cQ0E0T .* uNormalQ0E0T, problemData.N, problemData.basesOnQuad);
problemData.globKmuD = assembleVecEdgeMuFuncCont(problemData.g, problemData.g.markE0TbdrD, ...
                            problemData.cDCont, 2 * problemData.p + 1, problemData.basesOnGamma);

% Assemble outflow boundary conditions.
problemData.globSout = assembleMatEdgePhiIntMuVal(problemData.g, problemData.g.markE0TbdrN, ...
                                                  problemData.hatS, uNormalQ0E0T);

% %% HDG stuff
% %Evaluate u (=transport velocities) on every element
% problemData.uEval = evalFuncContAtEveryIntPoint(problemData.g, @(x1,x2) problemData.fluxCont( x1 ,x2, 1.), ...
%                                      problemData.N);
%                                  
% %Evaluate c (=solution) on every edge
% problemData.cEdge = evalFuncContAtEveryEdgeIntPoint( problemData.g, @(x1, x2) problemData.cDCont( x1 ,x2), ...
%                                          problemData.Nmu);
% % Evaluate advection velocity on every element
% problemData.uEdge = evalUContAtEveryEdgeIntPoint(problemData.g, @(x1, x2) problemData.fluxCont( x1 ,x2, 1.), ...
%                                     problemData.Nmu);
% %Evaluate f (= source term) on every element
% problemData.srcEval = evalSourceContAtEveryIntPoint(problemData.g, @(x1,x2) problemData.fCont( x1 ,x2 ), ...
%                                      problemData.N);
% problemData.globH = assembleVecElemPhiSource(problemData.g, problemData.N, ...
%                             problemData.srcEval, problemData.basesOnQuad);
% 
%                                                                                                 
% %Evaluate the flux on every edge
% problemData.fluxEdge = ...
%     evalFluxContAtEveryEdgeIntPoint( problemData.g, ...
%                                      problemData.g.markE0TbdrD, ...
%                                      @(x1, x2, c) problemData.fluxCont( x1 ,x2, c), ...
%                                      problemData.cEdge, problemData.Nmu);
%                         
% % Evaluate Dirichlet boundary condition for the first equation.
% problemData.globFphiD = assembleVecEdgePhiIntFlux( problemData.g, problemData.N, ...
%     problemData.fluxEdge, problemData.g.markE0TbdrD, problemData.basesOnQuad );
% 
% % 
% % problemData.globG = assembleMatElemPhiDphiFlux( problemData.g, problemData.N, problemData.uEval, problemData.hatG );
% u1Cont = @(x1, x2) problemData.u1Cont(x1, x2);
% u2Cont = @(x1, x2) problemData.u2Cont(x1, x2);
% problemData.globG = assembleMatElemDphiPhiFuncContVec(problemData.g, problemData.hatG, u1Cont, u2Cont);
% 
% % Flux on interior edges
% problemData.globS = assembleMatEdgePhiIntMuVal(problemData.g, problemData.g.markE0Tint, ...
%                                                  problemData.hatS, problemData.uEdge);
% % Outflow BC
% problemData.globSout = assembleMatEdgePhiIntMuVal(problemData.g, problemData.g.markE0TbdrN, ...
%                                                   problemData.hatS, problemData.uEdge);
% 
% % Assembly of Dirichlet boundary contributions
% % This has to be evaluated at t_new = t + dt!!
% problemData.globKmuD = assembleVecEdgeMuFuncContVal( problemData.g, problemData.g.markE0TbdrD, ...
%     @(x1,x2) problemData.cDCont( x1, x2), problemData.Nmu, problemData.basesOnGamma );
end % function
