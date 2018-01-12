% Performs all post-processing steps, such as error estimates of the final
% solution, etc.

%===============================================================================
%> @file darcyVert_sweVert/postprocessProblem.m
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%===============================================================================
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%>
%> This routine is called after the main loop.
%>
%> It can include things such as error estimates, an output operation of
%> the final solution, etc.
%>
%> @param  problemData  A struct with problem parameters and solution
%>                      vectors. @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description and precomputed fields.
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
function problemData = postprocessProblem(problemData)
problemData.darcyData = problemData.darcySteps.postprocessProblem(problemData.darcyData);

if problemData.isCouplingSWE
  % Coupling term for vertical velocity component
  K = problemData.sweData.g.numT;
  N = problemData.darcyData.N;
  
  markAreaE0T = problemData.sweData.g.markE0TbdrCoupling(:, 1) .* problemData.sweData.g.areaE0T(:, 1);
  
  % Upper edge (2) in Darcy problem is coupled to lower edge (1) in SWE 
  % problem:
  % Darcy values are evaluated on edge 2 and integrated over edge 1 in SWE
  % grid data.
  [Q, W] = quadRule1D(problemData.qOrd);
  [Q1, Q2] = gammaMapTetra(2, Q);
  
  % Evaluate K in quadrature points (2x2 cell array of K_PM x R arrays)
  KQ0E0T = cellfun(@(Kij) Kij(problemData.tEnd, problemData.darcyData.g.mapRef2Phy(1, Q1, Q2), problemData.darcyData.g.mapRef2Phy(2, Q1, Q2)), problemData.darcyData.KCont, 'UniformOutput', false);
  
  % Evaluate q1, q2 in quadrature points (K_PM x R arrays)
  q1Disc = reshape(problemData.darcyData.sysY(1 : K*N), N, K)';
  q2Disc = reshape(problemData.darcyData.sysY(K*N+1 : 2*K*N), N, K)';
  q1Q0E0T = q1Disc * problemData.darcyData.basesOnQuad.phi1D{problemData.qOrd}(:, :, 2)';
  q2Q0E0T = q2Disc * problemData.darcyData.basesOnQuad.phi1D{problemData.qOrd}(:, :, 2)';
  
  % Compute combined values (K_PM x R arrays)
  u1CouplingQ0E0T = KQ0E0T{1,1} .* q1Q0E0T + KQ0E0T{1,2} .* q2Q0E0T;
  u2CouplingQ0E0T = KQ0E0T{2,1} .* q1Q0E0T + KQ0E0T{2,2} .* q2Q0E0T;
  
  JuCoupling = markAreaE0T .* ((problemData.gCoupling.markE0TE0T{1} * u1CouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1)));
  problemData.sweData.globJuCoupling{1} = reshape((JuCoupling .* problemData.sweData.g.nuE0T(:, 1, 1)).', K*N, 1);
  problemData.sweData.globJuCoupling{2} = reshape((JuCoupling .* problemData.sweData.g.nuE0T(:, 1, 2)).', K*N, 1);
  problemData.sweData.globJwCoupling = reshape((markAreaE0T .* problemData.sweData.g.nuE0T(:, 1, 2) .* ((problemData.gCoupling.markE0TE0T{1} * u2CouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1)))).', K*N, 1);
 end % if

% Coupling term for vertical velocity component
% t = problemData.tEnd;
% K = problemData.darcyData.g.numT;
% N = problemData.darcyData.N;
% KDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.darcyData.g, @(x1,x2) c(t,x1,x2), problemData.darcyData.qOrd, ...
%                        problemData.darcyData.globM, problemData.darcyData.basesOnQuad), problemData.darcyData.KCont, 'UniformOutput', false);
%                      
% globRcouple1 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{1,1});
% globRcouple2 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{1,2});
% 
% problemData.sweData.globJuCoupling = { -(globRcouple1{1} *  problemData.darcyData.sysY(1 : K*N) + globRcouple2{1} * problemData.darcyData.sysY(K*N+1 : 2*K*N)), ...
%                                        -(globRcouple1{2} *  problemData.darcyData.sysY(1 : K*N) + globRcouple2{2} * problemData.darcyData.sysY(K*N+1 : 2*K*N)) };
%                                      
% globRcouple1 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{2,1});
% globRcouple2 = assembleMatEdgeTetraPhiIntPhiExtFuncDiscExtNu(problemData.gCoupling, problemData.sweData.g.markE0TbdrCoupling, ...
%                     problemData.sweData.hatRoffdiag, KDisc{2,2});
%                   
% problemData.sweData.globJwCoupling = globRcouple1{2} * problemData.darcyData.sysY(1 : K*N) + globRcouple2{2} * problemData.darcyData.sysY(K*N+1 : 2*K*N);

% u1DCont = @(x1,x2) problemData.sweData.u1DCont(t,x1,x2);
% globJu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.sweData.g, problemData.sweData.g.markE0TbdrCoupling, u1DCont, problemData.sweData.N, problemData.sweData.qOrd, problemData.sweData.basesOnQuad2D);
% problemData.sweData.globJuCoupling = globJu;
% 
% u2DCont = @(x1,x2) problemData.sweData.u2DCont(t,x1,x2);
% globJw = assembleVecEdgeTetraPhiIntFuncContNu(problemData.sweData.g, problemData.sweData.g.markE0TbdrCoupling, u2DCont, problemData.sweData.N, problemData.sweData.qOrd, problemData.sweData.basesOnQuad2D);
% problemData.sweData.globJwCoupling = globJw{2};

problemData.sweData = problemData.sweSteps.postprocessProblem(problemData.sweData);
problemData.error = [ problemData.sweData.error, problemData.darcyData.error ];
end

