% First step of the four-part algorithm in the main loop. Do-nothing
% function in the advection solver.

%===============================================================================
%> @file ./cahn-hilliard/configureProblem.m
%>
%> @brief First step of the four-part algorithm in the main loop. 
%>        Do-nothing function in the advection solver.
%===============================================================================
%>
%> @brief First step of the four-part algorithm in the main loop.
%>        Do-nothing function in the advection solver.
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
%> The Advection problem requires substepping due to the Runge-Kutta method
%> (see solveStep() and @ref RAWFK2016 for details). Thus, no terms can be
%> assembled here and this routine does nothing.
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
function problemData = preprocessStep(problemData, nStep) %#ok<INUSD>
% Here comes the assembly of the Double-Well Potential Terms!
K = problemData.K;
N = problemData.N;

Cdisc = reshape(problemData.sysY(1:K*N), N, K)';
problemData.globEold = assembleVecElemPhiFuncContPhi(problemData.g, problemData.basesOnQuad, Cdisc, problemData.ConcavePsi);

t = problemData.tau * nStep;
fDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont(t,x1,x2), 2 * problemData.p + 1, ...
                                 problemData.hatM, problemData.basesOnQuad);
problemData.globRHS = problemData.globM * reshape(fDisc', problemData.K * problemData.N, 1);

gNCont = @(x1,x2) problemData.gNContCon(t,x1,x2);
problemData.globKNCon = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                      gNCont, problemData.N, problemData.basesOnQuad, 2 * problemData.p + 1);

gNContPot = @(x1,x2) problemData.gNContPot(t,x1,x2);
problemData.globKNPot = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
                      gNContPot, problemData.N, problemData.basesOnQuad, 2 * problemData.p + 1);

if problemData.isFluxLim
  problemData.globMTimesCold = problemData.globM * problemData.sysY(1:K*N);
  problemData.globBarMTimesCold = problemData.globBarM * problemData.sysY(1:K*N);
end

end % function
