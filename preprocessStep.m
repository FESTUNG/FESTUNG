% First step of the four-part algorithm in the main loop.

%===============================================================================
%> @file darcy_swe_2dv/preprocessStep.m
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
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
t = (nStep-1) * problemData.darcyData.tau;

if problemData.isCouplingDarcy
  % Reset water height for coupling
  problemData.hSWE = zeros(size(problemData.sweData.cDiscRK{1, 1}));
  problemData.cSWE = cellfun(@(c) zeros(size(c)), problemData.sweData.cDiscRK(1, :), 'UniformOutput', false);
end % if

if problemData.isCouplingSWE
  % Coupling term for vertical velocity component
  K = problemData.sweData.g.numT;
  N = problemData.darcyData.N;
  
  markAreaE0T = problemData.sweData.g.markE0TbdrCoupling(:, 1) .* problemData.sweData.g.areaE0T(:, 1);
  
  % Upper edge (2) in Darcy problem is coupled to lower edge (1) in SWE problem:
  % Darcy values are evaluated on edge 2 and integrated over edge 1 in SWE grid data.
  [Q, W] = quadRule1D(problemData.qOrd);
  [Q1, Q2] = gammaMapTetra(2, Q);
  
  % Evaluate K in quadrature points (2x2 cell array of K_PM x R arrays)
  KQ0E0T = cellfun(@(Kij) Kij(t, problemData.darcyData.g.mapRef2Phy(1, Q1, Q2), problemData.darcyData.g.mapRef2Phy(2, Q1, Q2)), problemData.darcyData.KCont, 'UniformOutput', false);
  
  % Evaluate q1, q2 in quadrature points (K_PM x R arrays)
  q1Disc = reshape(problemData.darcyData.sysY(1 : K*N), N, K)';
  q2Disc = reshape(problemData.darcyData.sysY(K*N+1 : 2*K*N), N, K)';
  q1Q0E0T = q1Disc * problemData.darcyData.basesOnQuad.phi1D{problemData.qOrd}(:, :, 2).';
  q2Q0E0T = q2Disc * problemData.darcyData.basesOnQuad.phi1D{problemData.qOrd}(:, :, 2).';
  
  % Compute combined values (K_PM x R arrays)
  u1CouplingQ0E0T = KQ0E0T{1,1} .* q1Q0E0T + KQ0E0T{1,2} .* q2Q0E0T;
  u2CouplingQ0E0T = KQ0E0T{2,1} .* q1Q0E0T + KQ0E0T{2,2} .* q2Q0E0T;
  u1u1CouplingQ0E0T = u1CouplingQ0E0T .* u1CouplingQ0E0T;
  u1u2CouplingQ0E0T = u1CouplingQ0E0T .* u2CouplingQ0E0T;
  
  JuCoupling = markAreaE0T .* ((problemData.gCoupling.markE0TE0T{1} * u1CouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1)));
  problemData.sweData.globJuCoupling{1} = reshape((JuCoupling .* problemData.sweData.g.nuE0T(:, 1, 1)).', K*N, 1);
  problemData.sweData.globJuCoupling{2} = reshape((JuCoupling .* problemData.sweData.g.nuE0T(:, 1, 2)).', K*N, 1);
  problemData.sweData.globJwCoupling = reshape((markAreaE0T .* problemData.sweData.g.nuE0T(:, 1, 2) .* ((problemData.gCoupling.markE0TE0T{1} * u2CouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1)))).', K*N, 1);
  problemData.sweData.globJuuCoupling = reshape((markAreaE0T .* problemData.sweData.g.nuE0T(:, 1, 1) .* ((problemData.gCoupling.markE0TE0T{1} * u1u1CouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1)))).', K*N, 1);
  problemData.sweData.globJuwCoupling = reshape((markAreaE0T .* problemData.sweData.g.nuE0T(:, 1, 2) .* ((problemData.gCoupling.markE0TE0T{1} * u1u2CouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1)))).', K*N, 1);
end % if
end % function
