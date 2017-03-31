% First step of the four-part algorithm in the main loop.

%===============================================================================
%> @file darcyVert_sweVert/preprocessStep.m
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
t = nStep * problemData.darcyData.tau;

% Reset water height for coupling
problemData.hSWE = zeros(size(problemData.sweData.cDiscRK{1, 1}));

% Coupling term for vertical velocity component
K = problemData.darcyData.g.numT;
N = problemData.darcyData.N;
KDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.darcyData.g, @(x1,x2) c(t,x1,x2), problemData.darcyData.N, problemData.darcyData.qOrd, ...
                       problemData.darcyData.globM, problemData.darcyData.basesOnQuad), problemData.darcyData.KCont, 'UniformOutput', false);
                     
% globRcouple1 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{1,1});
% globRcouple2 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{1,2});
% 
% problemData.sweData.globJuCoupling = { globRcouple1{1} *  problemData.darcyData.sysY(1 : K*N) + globRcouple2{1} * problemData.darcyData.sysY(K*N+1 : 2*K*N), ...
%                                        globRcouple1{2} *  problemData.darcyData.sysY(1 : K*N) + globRcouple2{2} * problemData.darcyData.sysY(K*N+1 : 2*K*N) };
                                     
% globRcouple1 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{2,1});
% globRcouple2 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{2,2});
%                   
% problemData.sweData.globJwCoupling = globRcouple1{2} * problemData.darcyData.sysY(1 : K*N) + globRcouple2{2} * problemData.darcyData.sysY(K*N+1 : 2*K*N);

% globRcouple1 = assembleMatEdgeTetraPhiIntPhiIntFuncDiscExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, problemData.hatS, KDisc{1,1}, reshape(problemData.darcyData.sysY(1 : K*N), N, K).');
% globRcouple2 = assembleMatEdgeTetraPhiIntPhiIntFuncDiscExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, problemData.hatS, KDisc{1,2}, reshape(problemData.darcyData.sysY(K*N+1 : 2*K*N), N, K).');
end % function
