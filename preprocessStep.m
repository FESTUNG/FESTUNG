% First step of the four-part algorithm in the time stepping loop.

%===============================================================================
%> @file
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
%>  1. darcy_2dv/preprocessStep.m
%>  2. darcy_2dv/solveStep.m
%>  3. darcy_2dv/postprocessStep.m
%>  4. darcy_2dv/outputStep.m
%> 
%> This routine is executed first in each loop iteration.
%> It takes care of the assembly of time-dependent matrices and right hand
%> side vectors. Furthermore, it evaluates the coefficient functions and
%> boundary conditions for the current time step.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      darcy_2dv/configureProblem.m and 
%>                      darcy_2dv/preprocessProblem.m. @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct without any modifications.
%>                      @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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

%% L2-projections of algebraic coefficients.
if iscell(problemData.DCont)
  DDisc = cellfun(@(c) projectFuncCont2DataDiscQuadri(problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
            problemData.globM, problemData.basesOnQuad), problemData.DCont, 'UniformOutput', false);
else
  DDisc = projectFuncCont2DataDiscQuadri(problemData.g, @(x1,x2) problemData.DCont(t,x1,x2), problemData.qOrd, ...
            problemData.globM, problemData.basesOnQuad);
end % if
fDisc = projectFuncCont2DataDiscQuadri(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), problemData.qOrd, ...
          problemData.globM, problemData.basesOnQuad);
                             
%% Assembly of time-dependent global matrices.
problemData.globG = assembleMatElemDphiPhiFuncDisc(problemData.g, problemData.hatG, DDisc);
problemData.globR = assembleMatEdgePhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, ...
                      problemData.hatRdiag, problemData.hatRoffdiag, DDisc);

%% Assembly of Dirichlet boundary contributions.
hDCont = @(x1,x2) problemData.hDCont(t,x1,x2);
problemData.globRD = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(problemData.g, ...
                      problemData.g.markE0TbdrD | problemData.g.markE0TbdrCoupling, problemData.hatRdiag, DDisc);
problemData.globJD = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrD, ...
                      hDCont, problemData.N, problemData.basesOnQuad, problemData.qOrd);
problemData.globKD = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrD, ...
                      hDCont, problemData.N, problemData.basesOnQuad, problemData.qOrd, ones(problemData.g.numT, 4));
                  
%% Assembly of Neumann boundary contributions.
gNCont = @(x1,x2) problemData.gNCont(t,x1,x2);
problemData.globKN = assembleVecEdgePhiIntFuncCont(problemData.g, ...
                     problemData.g.markE0TbdrN, gNCont, problemData.N, problemData.basesOnQuad, problemData.qOrd);
                   
%% Assembly of the source contribution.
problemData.globL = problemData.globM * reshape(fDisc', problemData.g.numT * problemData.N, 1);

end % function
