% Last step of the three-part algorithm in the iterateSubSteps loop of each step
% in the main loop.

%===============================================================================
%> @file
%>
%> @brief Last step of the three-part algorithm in the iterateSubSteps loop of 
%>        each step in the main loop.
%===============================================================================
%>
%> @brief Last step of the three-part algorithm in the iterateSubSteps loop of 
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
%> This routine is executed last in each substep loop iteration and is intended
%> to post-process the solution computed by solveSubStep().
%>
%> @param  pd           A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval pd           The input struct enriched with post-processed data
%>                      the next substep. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>                      Modified 08/23/16 by Hennes Hajduk
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
function pd = postprocessSubStep(pd, nStep, nSubStep)
K = pd.K;
N = pd.N;

% Reshape linearized vector to solution vectors
pd.cDisc(:,:,1) = reshape(pd.cDiscRK(        1 :   K*N), N, K).';
pd.cDisc(:,:,2) = reshape(pd.cDiscRK(  K*N + 1 : 2*K*N), N, K).';
pd.cDisc(:,:,3) = reshape(pd.cDiscRK(2*K*N + 1 : 3*K*N), N, K).';

% Slope Limiting
if pd.isSlopeLim
  pd = applySlopeLimiter(pd, pd.slopeLimName, nSubStep);
end % if

% Ensure water height doesn't fall below threshold
pd.cDisc(:,:,1) = correctMinValueExceedanceDisc(pd.cDisc(:,:,1), pd.sysMinValueCorrection, nStep, pd.zbV0T + pd.minTol, pd.elevTol);
pd.cDiscRK(1 : pd.K*pd.N) = reshape(pd.cDisc(:,:,1).', N*K, 1);

pd.isSubSteppingFinished = nSubStep >= length(pd.tLvls);
end % function
