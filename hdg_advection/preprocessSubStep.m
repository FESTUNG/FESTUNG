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

problemData.timeRK = problemData.t+ problemData.tabRK.C(nSubStep) * problemData.dt;

problemData.cDiscRkStep = zeros( K * N, 1 );
% problemData.cDiscRkStep = reshape( problemData.cDisc', size(problemData.globM, 1), 1 );
% for i=1:nSubStep-1
%     problemData.cDiscRkStep = problemData.cDiscRkStep + problemData.dt * problemData.tabRK.A(nSubStep,i) .* problemData.cDiscRK{i};
% end

problemData.cDiscRkRHS = zeros( K * N, 1 );
% RK RHS
for i=1:nSubStep-1
    problemData.cDiscRkRHS = problemData.cDiscRkRHS + problemData.tabRK.A(nSubStep,i) .* problemData.cDiscRK{i};
end

%% HDG stuff
%Evaluate u (=transport velocities) on every element
problemData.uEval(:,:,1) = evalFuncContAtEveryIntPoint(problemData.g, @(x1,x2) problemData.u1Cont( problemData.timeRK, x1,x2), ...
                                     problemData.N, problemData.basesOnQuad);
problemData.uEval(:,:,2) = evalFuncContAtEveryIntPoint(problemData.g, @(x1,x2) problemData.u2Cont( problemData.timeRK,x1,x2), ...
                                     problemData.N, problemData.basesOnQuad);

%Evaluate c (=solution) on every edge
problemData.cEdge = evalFuncContAtEveryEdgeIntPoint( problemData.g, @(x1, x2) problemData.cDCont(  problemData.timeRK, x1 ,x2), ...
                                         problemData.Nlambda);
                                     
%Evaluate the flux on every edge
problemData.fluxEdge = evalFluxContAtEveryEdgeIntPoint(problemData.g, problemData, @(x1, x2, c) problemData.fluxCont( problemData.timeRK, x1 ,x2, c), ...
                                           problemData.cEdge, problemData.Nlambda);
                               
% Term III.4
problemData.globFgamma = assembleVecEdgePhiIntFlux( problemData.g, problemData.N, ...
    problemData.fluxEdge, problemData.g.markE0TbdrD, problemData.basesOnQuad );
% Term III.6
problemData.globCd = assembleVecEdgePhiIntVal( problemData.g, problemData.N, ...
                                               problemData.cEdge, problemData.g.markE0TbdrD, ...
                                               problemData.basesOnQuad );

% Term II
problemData.globG = assembleMatElemPhiDphiFlux( problemData.g, problemData.N, problemData.uEval, problemData.hatGbarOnQuad );
               
% Evaluate advection velocity on every element
problemData.uEdge = evalUContAtEveryEdgeIntPoint(problemData.g, @(x1, x2, c) problemData.fluxCont( problemData.timeRK, x1 ,x2, 1.), ...
                                    problemData.Nlambda);

problemData.globS = assembleMatEdgeMuPhiIntFlux( problemData.g, problemData.g.markE0Tint, ...
                                                 problemData.uEdge, problemData.hatSbarOnQuad );
% Neumann BC
problemData.globSN = assembleMatEdgeMuPhiIntFlux( problemData.g, problemData.g.markE0TbdrN, ...
                                                  problemData.uEdge, problemData.hatSbarOnQuad );
% Assembly of Dirichlet boundary contributions
% This has to be evaluated at t_new = t + dt!!
problemData.globKDlambda = assembleVecEdgeMuFuncContVal( problemData.g, problemData.g.markE0TbdrD, ...
    @(x1,x2) problemData.cDCont( problemData.timeRK, x1, x2), problemData.Nlambda, problemData.basesOnGamma );

% Reshape cDisc to have a vector
% problemData.cDiscReshaped = reshape( problemData.cDisc', size(problemData.globM, 1), 1 );
% problemData.cDiscLambdaReshaped = reshape( problemData.cDiscLambda', size(problemData.globP, 1), 1 );
% M*cDisc, should I store it?

end % function
