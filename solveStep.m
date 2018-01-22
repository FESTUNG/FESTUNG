% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file darcy_swe_2dv/solveStep.m
%>
%> @brief Second step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
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
%> This routine is executed second in each loop iteration and is intended to
%> produce the solution at the next step, e.g., at a new time-level.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next step. @f$[\text{struct}]@f$
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
function problemData = solveStep(problemData, nStep)
%% Solve SWE time steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep);

%% Solve Darcy time step
problemData.darcyData = problemData.darcySteps.preprocessStep(problemData.darcyData, nStep);

if problemData.isCouplingDarcy
%   % Compute coupling terms
%   problemData.hSWE = reshape(1 / problemData.darcyData.tau * problemData.hSWE.', [], 1);  % mass conservative
% %   problemData.hSWE = reshape(problemData.sweData.cDiscRK{end, 1}.', [], 1);  % analytically correct
%   for m = 1 : 2
%     problemData.darcyData.globJcouple{m} = problemData.darcyData.tildeGlobQcouple{m} * problemData.hSWE;
%   end % for
%   problemData.darcyData.globKcouple = problemData.darcyData.tildeGlobScouple * problemData.hSWE;
  
  % Coupling term for Darcy head
  K = problemData.darcyData.g.numT;
  N = problemData.sweData.N;
  
  % Upper edge (2) in Darcy problem is coupled to lower edge (1) in SWE problem:
  % SWE values are evaluated on edge 1 and integrated over edge 2 in Darcy grid data.
  [~, W] = quadRule1D(problemData.qOrd);
  
  % Evaluate primary variables in quadrature points of bottom edge of SWE domain
  hQ0E0T = (1/problemData.darcyData.tau) * problemData.cSWE{1} * problemData.sweData.basesOnQuad1D.phi1D{problemData.qOrd}.';
  u1Q0E0T = (1/problemData.darcyData.tau) * problemData.cSWE{2} * problemData.sweData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 2).';
  hCouplingQ0E0T = problemData.markT2DT * hQ0E0T + 0.5 * ( u1Q0E0T .* u1Q0E0T );
  
  % Integrate boundary condition
  markAreaNuE0T = bsxfun(@times, problemData.darcyData.g.markE0TbdrCoupling(:, 2) .* problemData.darcyData.g.areaE0T(:, 2), squeeze(problemData.darcyData.g.nuE0T(:, 2, :)));
  globInt = (problemData.markE0TE0T.' * hCouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.darcyData.basesOnQuad.phi1D{problemData.qOrd}(:, :, 2));
  
  %globJcouple = cell(2,1);
  % Apply coefficients
  for m = 1 : 2
    problemData.darcyData.globJcouple{m} = reshape((markAreaNuE0T(:, m) .* globInt).', K*N, 1);
  end % for m
  problemData.darcyData.globKcouple = reshape((problemData.darcyData.g.markE0TbdrCoupling(:, 2) .* globInt).', K*N, 1);
end % if
                   
problemData.darcyData = problemData.darcySteps.solveStep(problemData.darcyData, nStep);
problemData.darcyData = problemData.darcySteps.postprocessStep(problemData.darcyData, nStep);
end % function
