% First step of the three-part algorithm in the iterateSubSteps loop of each
% step in the main loop.

%===============================================================================
%> @file transport/preprocessSubStep.m
%>
%> @brief First step of the three-part algorithm in the iterateSubSteps loop of 
%>        each step in the main loop.
%===============================================================================
%>
%> @brief First step of the three-part algorithm in the iterateSubSteps loop of 
%>        each step in the main loop.
%>
%> The iterateSubSteps loop repeatedly executes three steps until the number of
%> substep iterations equals the order of the underlying Runge-Kutta method.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%> 
%> This routine is executed first in each substep loop iteration and is intended
%> to execute preprocessing operations, e.g., evaluate boundary conditions or
%> right hand side values, assemble time-dependent matrix blocks, etc.
%>
%> @param  pd           A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval pd           The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = preprocessSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;
p = problemData.p;
qOrd2D = max(2*p,1);

% Project velocity to DG space, if not given
if ~problemData.isCoupling
  problemData.uHDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.uHCont(problemData.timeLvls(nSubStep),x1,x2), 2*p, ...
                                                problemData.hatM, problemData.basesOnQuad);
  problemData.vHDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.vHCont(problemData.timeLvls(nSubStep),x1,x2), 2*p, ...
                                                problemData.hatM, problemData.basesOnQuad);
  % Evaluate normal velocity in quadrature points of edges
  problemData.vNormalOnQuadEdge = computeFuncContNuOnQuadEdge(problemData.g, @(x1,x2) problemData.uHCont(problemData.timeLvls(nSubStep),x1,x2), ...
                                                              @(x1,x2) problemData.vHCont(problemData.timeLvls(nSubStep),x1,x2), 2*p+1);
end % if

% Evaluate cDisc in all quadrature points
problemData.cHQ0T = cellfun(@(c) (reshape(c, N, K).' * problemData.basesOnQuad.phi2D{qOrd2D}.'), problemData.cDiscRK, 'UniformOutput', false);

end % function
