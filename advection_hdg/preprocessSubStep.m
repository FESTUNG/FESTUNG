% Preprocessing of the Runge-Kutta step.

%===============================================================================
%> @file
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
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
%> @author Alexander Jaust, 2017.
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
function problemData = preprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
% Select time level for current Runge-Kutta stage
t = problemData.t(nSubStep);
cDCont = @(x1, x2) problemData.cDCont(t, x1 ,x2);
u1Cont = @(x1, x2) problemData.u1Cont(t, x1, x2);
u2Cont = @(x1, x2) problemData.u2Cont(t, x1, x2);
fCont = @(x1,x2) problemData.fCont(t, x1, x2);

% Evaluate normal advection velocity on every edge.
uNormalQ0E0T = computeFuncContNuOnQuadEdge(problemData.g, u1Cont, u2Cont, problemData.qOrd);

% Determine inflow- and outflow edges
[~, W] = quadRule1D(problemData.qOrd);
markE0TbdrOut = problemData.g.markE0Tbdr & sum(bsxfun(@times, reshape(W, 1, 1, []), uNormalQ0E0T), 3) > 0;
markE0TbdrIn = problemData.g.markE0Tbdr & ~markE0TbdrOut;

% Assemble source term.
problemData.globH = problemData.globMphi * reshape(projectFuncCont2DataDisc(problemData.g, fCont, ...
                        problemData.qOrd, problemData.hatM, problemData.basesOnQuad)', [], 1);

% Assemble element integral contributions.
problemData.globG = assembleMatElemDphiPhiFuncContVec(problemData.g, problemData.hatG, ...
                        u1Cont, u2Cont, problemData.qOrd);

% Assemble flux on interior and outflow edges
problemData.globS = assembleMatEdgePhiIntMuVal(problemData.g, problemData.g.markE0Tint | markE0TbdrOut, ...
                        problemData.hatS, uNormalQ0E0T);
problemData.globKmuOut = assembleMatEdgePhiIntMu(problemData.g, markE0TbdrOut, problemData.hatRmu)';
                      
% Assemble inflow boundary conditions.
cQ0E0T = computeFuncContOnQuadEdge(problemData.g, cDCont, problemData.qOrd);
problemData.globFphiIn = assembleVecEdgePhiIntVal(problemData.g, markE0TbdrIn, ...
                            cQ0E0T .* uNormalQ0E0T, problemData.N, problemData.basesOnQuad);
problemData.globKmuIn = assembleVecEdgeMuFuncCont(problemData.g, markE0TbdrIn, ...
                            cDCont, problemData.basesOnQuad, problemData.qOrd);
end % function
