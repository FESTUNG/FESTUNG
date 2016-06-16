% Third step of the four-part algorithm in the main loop.

%===============================================================================
%> @file template/postprocessStep.m
%>
%> @brief Third step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Third step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the number of
%> iterations provided by configureProblem in the parameter
%> <code>numSteps</code> is reached. These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed third in each loop iteration and is intended to
%> post-process the solution computed by solveStep().
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures, as provided 
%>                      by configureProblem() and preprocessProblem(). 
%>                      @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with post-processed data
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
function pd = postprocessStep(pd, nStep)
% Ensure water height doesn't fall below threshold
pd.cDisc(:,:,1) = correctMinValueExceedanceDisc(pd.cDisc(:,:,1), pd.sysMinValueCorrection, nStep, pd.zbLagr + pd.minTol, 20);

% Update time level and check for simulation end
pd.t = pd.t + pd.dt;
pd.isFinished = pd.t >= pd.tEnd;

% Check for steady state convergence
if pd.isSteadyState && pd.changeL2 < pd.convergenceCriterion
  fprintf('Steady state is reached.\n')
  pd.isFinished = true;
end % if
end

