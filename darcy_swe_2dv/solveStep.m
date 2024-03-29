% Second step of the four-part algorithm in the main loop. It solves a
% number of free-flow time steps and one subsurface time step.

%===============================================================================
%> @file
%>
%> @brief Second step of the four-part algorithm in the main loop. It solves a
%>        number of free-flow time steps and one subsurface time step.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop. It solves a
%>        number of free-flow time steps and one subsurface time step.
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. darcy_swe_2dv/preprocessStep.m
%>  2. darcy_swe_2dv/solveStep.m
%>  3. darcy_swe_2dv/postprocessStep.m
%>  4. darcy_swe_2dv/outputStep.m
%> 
%> This routine carries out substepping to solve a number of free flow time
%> steps and afterwards computes the flux from free flow domain to
%> subsurface domain before executing one time step of the subsurface
%> problem.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by 
%>                      darcy_swe_2dv/configureProblem.m  and 
%>                      darcy_swe_2dv/preprocessProblem.m @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
function problemData = solveStep(problemData, nStep)
%% Solve SWE time steps
problemData.isSubSteppingFinished = false;
problemData = iterateSubSteps(problemData, nStep, problemData.subStepHandles);

%% Solve Darcy time step
problemData.darcyData = problemData.darcySteps.preprocessStep(problemData.darcyData, nStep);

% Coupling term for Darcy head
if problemData.isCouplingDarcy
  qOrd = problemData.qOrd;
  
  K = problemData.darcyData.g.numT;
  N = problemData.darcyData.N;
  
  % Upper edge (2) in Darcy problem is coupled to lower edge (1) in SWE problem:
  % SWE values are evaluated on edge 1 and integrated over edge 2 in Darcy grid data.
  [~, W] = quadRule1D(qOrd);
  
  % Integrate boundary condition
  markAreaNuE0T = bsxfun(@times, problemData.darcyData.g.markE0TbdrCoupling(:, 2) .* problemData.darcyData.g.areaE0T(:, 2), squeeze(problemData.darcyData.g.nuE0T(:, 2, :)));
  globInt = (problemData.markE0TE0T.' * problemData.hCouplingQ0E0T) * (repmat(W(:), 1, N) .* problemData.darcyData.basesOnQuad.phi1D{qOrd}(:, :, 2));
  
  % Apply coefficients
  for m = 1 : 2
    problemData.darcyData.globJcouple{m} = reshape(bsxfun(@times, markAreaNuE0T(:, m), globInt).', K*N, 1);
  end % for m
  
  if problemData.darcyData.isJumpCoupling
    problemData.darcyData.globKcouple = reshape(bsxfun(@times, problemData.darcyData.g.markE0TbdrCoupling(:, 2), globInt).', K*N, 1);
  end % if
end % if
                   
problemData.darcyData = problemData.darcySteps.solveStep(problemData.darcyData, nStep);
problemData.darcyData = problemData.darcySteps.postprocessStep(problemData.darcyData, nStep);
end % function
