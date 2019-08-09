% First step of the four-part algorithm in the main loop.

%===============================================================================
%> @file
%>
%> @brief First step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief First step of the four-part algorithm in the main loop.
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
%> This routine is executed first in each loop iteration and is intended to
%> execute preprocessing operations, e.g., evaluate boundary conditions or
%> right hand side values, assemble time-dependent matrix blocks, etc.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
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
function problemData = preprocessStep(problemData, nStep)
%% Compute current time level and create function handles
t = problemData.t0 + problemData.dt * nStep;
dCont = @(x1,x2) problemData.dCont(t,x1,x2);
rCont = @(x1,x2) problemData.rCont(t,x1,x2);
v1Cont = @(x1,x2) problemData.v1Cont(t,x1,x2);
v2Cont = @(x1,x2) problemData.v2Cont(t,x1,x2);
uDCont = @(x1,x2) problemData.uDCont(t,x1,x2);
gNCont = @(x1,x2) problemData.gNCont(t,x1,x2);
fCont = @(x1,x2) problemData.fCont(t,x1,x2);
%% Assemble advection and reaction matrices/vectors
% Evaluate normal velocity in quadrature points of edges
vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, v1Cont, v2Cont, 2*problemData.p+1);
problemData.globAadv = assembleMatElemDphiPhiFuncContVec(problemData.g, problemData.refElemDphiPhiPerQuad, v1Cont, v2Cont);
problemData.globAreac = assembleMatElemPhiPhiFuncCont(problemData.g, problemData.refElemPhiPhiPerQuad, rCont);
problemData.globBadv = assembleMatEdgePhiPhiValUpwind(problemData.g, true(problemData.g.numT, 3), ...
  problemData.refEdgePhiIntPhiIntPerQuad, problemData.refEdgePhiIntPhiExtPerQuad, vNormalOnQuadEdge);
problemData.globJadv = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
  uDCont, vNormalOnQuadEdge, problemData.N, problemData.basesOnQuad);
%% Assemble right hand side vector
problemData.sysF = reshape(projectFuncCont2DataDisc(problemData.g, fCont, 2*problemData.p, ...
                           problemData.refElemPhiPhi, problemData.basesOnQuad)', [], 1);
%% Assemble diffusion matrices/vectors
problemData.globJN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, gNCont, ...
  problemData.N, problemData.basesOnQuad);
problemData.globJjmp = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrD, uDCont, ...
  problemData.N, problemData.basesOnQuad, 2 * problemData.p + 1, ones(problemData.g.numT, 3));
if problemData.isIP  % IP-specific matrices/vectors
  problemData.globAIP = assembleMatElemDphiDphiFuncCont(problemData.g, problemData.refElemDphiDphiPerQuad, dCont);
  problemData.globAIP = problemData.globAIP{1} + problemData.globAIP{2};
  problemData.globBIP = assembleMatEdgeDphiPhiFuncContNu(problemData.g, problemData.g.markE0Tint, ...
    problemData.refEdgeDphiIntPhiIntPerQuad, problemData.refEdgeDphiIntPhiExtPerQuad, dCont);
  problemData.globBIP = 0.5 * (problemData.globBIP{1}' + problemData.globBIP{2}');
  problemData.globBIPD = assembleMatEdgeDphiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                                                          problemData.refEdgeDphiIntPhiIntPerQuad, dCont);
  problemData.globBIPD = 0.5 * (problemData.globBIPD{1}' + problemData.globBIPD{2}');
  problemData.globBsym = assembleMatEdgeDphiPhiFuncContNu(problemData.g, ~problemData.g.markE0TbdrN, ...
    problemData.refEdgeDphiIntPhiIntPerQuad, problemData.refEdgeDphiIntPhiExtPerQuad, dCont);
  problemData.globBsym = problemData.globBsym{1} + problemData.globBsym{2};
  problemData.globJsym = assembleVecEdgeDphiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
    @(x1,x2) dCont(x1,x2) .* uDCont(x1,x2), problemData.N, problemData.basesOnQuad);
  problemData.globJsym = problemData.globJsym{1} + problemData.globJsym{2};
else                 % LDG-specific matrices/vectors
  problemData.globAu = assembleMatElemDphiPhiFuncContVec(problemData.g, problemData.refElemDphiPhiPerQuad, dCont, dCont);
  problemData.globBu = assembleMatEdgePhiPhiFuncContNu(problemData.g, problemData.g.markE0Tint, ...
    problemData.refEdgePhiIntPhiIntPerQuad, problemData.refEdgePhiIntPhiExtPerQuad, dCont);
  problemData.globBuD = assembleMatEdgePhiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                                                              problemData.refEdgePhiIntPhiIntPerQuad, dCont);
  problemData.globJD = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, uDCont, ...
                                                       problemData.N, problemData.basesOnQuad);
end
end % function
