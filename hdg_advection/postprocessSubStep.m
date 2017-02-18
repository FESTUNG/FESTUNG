% Third step of the three-part substepping algorithm for each Runge-Kutta stage.

%===============================================================================
%> @file advection/postprocessSubStep.m
%>
%> @brief Third step of the three-part substepping algorithm for each Runge-Kutta stage.
%===============================================================================
%>
%> @brief Third step of the three-part substepping algorithm for each Runge-Kutta stage.
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
%> This routine is executed third in each loop iteration.
%> It decides whether the substepping is finished and updates
%> <code>problemData.isSubSteppingFinished</code> accordingly.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop.
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with postprocessed data
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
function problemData = postprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
% problemData.isSubSteppingFinished = nSubStep >= length(problemData.omega);

stab = problemData.stab;
diagRK = 1.;


problemData.cDiscReshaped = reshape( problemData.cDisc', size(problemData.globM, 1), 1 );
problemData.cDiscLambdaReshaped = reshape( problemData.cDiscLambda', size(problemData.globP, 1), 1 );
%Evaluate the flux on every edge
problemData.fluxEdge = ...
    evalFluxContAtEveryEdgeIntPoint(problemData.g, problemData, ...
                         @(x1, x2, c) problemData.fluxCont( problemData.timeRK, x1 ,x2, c), ...
                                           problemData.cEdge, problemData.Nlambda);
                               
% Term III.4
problemData.globFgamma = assembleVecEdgePhiIntFlux( problemData.g, problemData.N, problemData.fluxEdge, problemData.g.markE0TbdrD, problemData.basesOnQuad );
% Term III.6
problemData.globCd = assembleVecEdgePhiIntVal( problemData.g, problemData.N, ...
                                               problemData.cEdge, problemData.g.markE0TbdrD, ...
                                               problemData.basesOnQuad );
% Term II
problemData.globG = assembleMatElemPhiDphiFlux( problemData.g, problemData.N, problemData.uEval, problemData.hatGbarOnQuad );

problemData.cDiscRK{nSubStep} = - diagRK .* stab * problemData.globFgamma - diagRK .* stab * problemData.globCd ...
                                - (- diagRK .* problemData.globG{1} - diagRK .* problemData.globG{2} ...
                                + diagRK .* stab * problemData.globRphi) * problemData.cDiscReshaped ...
                                - (problemData.globS{1} + problemData.globS{2} - stab * problemData.globRlambda + problemData.globSN{1} + problemData.globSN{2}  + stab * problemData.globRD) * problemData.cDiscLambdaReshaped;

                            
if (nSubStep == problemData.tabRK.s)
    problemData.isSubSteppingFinished = true;
end
end % function
