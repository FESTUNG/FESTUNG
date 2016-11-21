% Last step of the three-part algorithm in the iterateSubSteps loop of each step
% in the main loop.

%===============================================================================
%> @file template/postprocessSubStep.m
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
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval problemData  The input struct enriched with post-processed data
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
function problemData = postprocessSubStep(problemData, nStep, nSubStep)
K = problemData.K;
N = problemData.N;

% Reshape linearized vector to solution vectors
problemData.cDisc(:,:,1) = reshape(problemData.cDiscRK(        1 :   K*N), N, K).';
problemData.cDisc(:,:,2) = reshape(problemData.cDiscRK(  K*N + 1 : 2*K*N), N, K).';
problemData.cDisc(:,:,3) = reshape(problemData.cDiscRK(2*K*N + 1 : 3*K*N), N, K).';

for i = 1 : length(problemData.slopeLimList)
  switch problemData.slopeLimList{i}
    case 'xi'
      problemData.cDisc(:,:,1) = applySlopeLimiterDisc(problemData.g, problemData.cDisc(:,:,1), problemData.g.markV0TbdrD, problemData.ramp(problemData.tLvls(nSubStep)/86400) * ...
                                  (problemData.dataV0Triv + problemData.dataV0Tos), problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim); % TODO time
    otherwise
      error('Slope limiting not implemented for variables other than free surface elevation.')
  end % switch
end % for

% Ensure water height doesn't fall below threshold
problemData.cDisc(:,:,1) = correctMinValueExceedanceDisc(problemData.cDisc(:,:,1), problemData.sysMinValueCorrection, nStep, problemData.zbLagr + problemData.minTol, problemData.elevTol);
problemData.cDiscRK(1 : problemData.K*problemData.N) = reshape(problemData.cDisc(:,:,1).', problemData.N * problemData.K, 1);

problemData.isSubSteppingFinished = nSubStep >= length(problemData.tLvls);
end % function
