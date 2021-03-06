% First step of the four-part algorithm in the time stepping loop.

%===============================================================================
%> @file advection-diffusion/preprocessStep.m
%>
%> @brief First step of the four-part algorithm in the time stepping loop. 
%===============================================================================
%>
%> @brief First step of the four-part algorithm in the time stepping loop. 
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
%> This routine is executed first in each loop iteration.
%> It takes care of the assembly of time-dependent matrices and right hand
%> side vectors. Furthermore, it evaluates the coefficient functions and
%> boundary conditions for the current time step.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct without any modifications.
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
function problemData = preprocessStep(problemData, nStep)
t = nStep * problemData.tau;
problemData.t = t;

K = problemData.K;
N = problemData.N;

%% L2-projections of algebraic coefficients.
dDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.dCont(t,x1,x2), 2 * problemData.p, ...
                                 problemData.hatM, problemData.basesOnQuad);
fDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), 2 * problemData.p, ...
                                 problemData.hatM, problemData.basesOnQuad);

u1Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u1Cont(t,x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);
u2Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.u2Cont(t,x1,x2), ...
                                  2*problemData.p, problemData.hatM, problemData.basesOnQuad);

%% Assembly of time-dependent global matrices.
problemData.globG = assembleMatElemDphiPhiFuncDisc(problemData.g, problemData.hatG, dDisc);
problemData.globR = assembleMatEdgePhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatRdiag, problemData.hatRoffdiag, dDisc, problemData.g.areaNuE0Tint);

% Evaluate normal velocity in quadrature points of edges
vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.u1Cont(t,x1,x2), ...
                      @(x1,x2) problemData.u2Cont(t,x1,x2), 2*problemData.p+1);

% Assembly of time-dependent global matrices
problemData.globGadv = assembleMatElemDphiPhiFuncDiscVec(problemData.g, problemData.hatG, u1Disc, u2Disc);
problemData.globRadv = assembleMatEdgePhiPhiValUpwind(problemData.g, ~problemData.g.markE0TbdrN, ...
                                                   problemData.hatRdiagOnQuad, problemData.hatRoffdiagOnQuad, ...
                                                   vNormalOnQuadEdge, problemData.g.areaE0TbdrNotN);

%% Assembly of Dirichlet boundary contributions.
cDCont = @(x1,x2) problemData.cDCont(t,x1,x2);
problemData.globRD = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrD, ...
                      problemData.hatRdiag, dDisc, problemData.g.areaNuE0TbdrD);
problemData.globJD = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                      cDCont, problemData.N, problemData.basesOnQuad, 2 * problemData.p + 1, problemData.g.areaNuE0TbdrD);
problemData.globKD = problemData.eta * assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrD, ...
                      cDCont, problemData.N, problemData.basesOnQuad, 2 * problemData.p + 1);
%% Assembly of Neumann boundary contributions.
problemData.globKN = assembleVecEdgePhiIntFuncDiscIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                      dDisc, @(x1,x2) problemData.gNCont(t,x1,x2), problemData.basesOnQuad, problemData.g.areaE0TbdrN);

% Assembly of Neumann boundary contributions
gNUpwind = @(x1,x2) (problemData.gNCont(t,x1,x2) <= 0) .* problemData.gNCont(t,x1,x2);
problemData.globKNadv = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                                                   gNUpwind, N, problemData.basesOnQuad);

%% Assembly of the source contribution.
problemData.globL = problemData.globM * reshape(fDisc', problemData.K * problemData.N, 1);
end % function