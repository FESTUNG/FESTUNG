% Last step of the three-part algorithm in the iterateSubSteps loop of each step
% in the main loop.

%===============================================================================
%> @file sweInverse/postprocessSubStep.m
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
function problemData = postprocessSubStep(problemData, nStep, nSubStep)
K = problemData.K;
N = problemData.N;

% Reshape linearized vector to solution vectors
problemData.cDisc(:,:,1) = reshape(problemData.cDiscRK(        1 :   K*N), N, K).';
problemData.cDisc(:,:,2) = reshape(problemData.cDiscRK(  K*N + 1 : 2*K*N), N, K).';
problemData.cDisc(:,:,3) = reshape(problemData.cDiscRK(2*K*N + 1 : 3*K*N), N, K).';

for i = 1 : length(problemData.slopeLimList)
  switch problemData.slopeLimList{i}
    case 'elevation'
      problemData.cDisc(:,:,1) = applySlopeLimiterDisc(problemData.g, problemData.cDisc(:,:,1), problemData.g.markV0TbdrD, ...
                                  problemData.ramp(getdefault(problemData.tLvls, nSubStep+1, problemData.t + problemData.dt)/86400) * (problemData.xiV0Triv + problemData.xiV0Tos), ...
                                  problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim);
    case 'momentum'
      problemData.cDisc(:,:,2) = applySlopeLimiterDisc(problemData.g, problemData.cDisc(:,:,2), problemData.g.markV0TbdrRI, ...
                                  problemData.ramp(getdefault(problemData.tLvls, nSubStep+1, problemData.t + problemData.dt)/86400) * problemData.uHV0Triv, ...
                                  problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim);
      
      problemData.cDisc(:,:,3) = applySlopeLimiterDisc(problemData.g, problemData.cDisc(:,:,3), problemData.g.markV0TbdrRI, ...
                                  problemData.ramp(getdefault(problemData.tLvls, nSubStep+1, problemData.t + problemData.dt)/86400) * problemData.vHV0Triv, ...
                                  problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim);
    otherwise
      error('Slope limiting not implemented for non primary variables.')
  end % switch
end % for

if ~problemData.isSteadyState
	hotstartData = readHotstart([problemData.elevationInput '_' num2str(nStep) '.mat']);
	assert(isstruct(hotstartData) && isfield(hotstartData, 'cDisc') && isfield (hotstartData, 't'), 'Hotstart data does not contain cDisc or t');
	problemData.xiDisc = addNoise(hotstartData.cDisc, problemData.noiseLvl);
	validateattributes(problemData.xiDisc, {'numeric'}, {'size', [K N]}, mfilename, 'problemData.xiDisc');
end % if

% Compute change
zbDisc = problemData.rampInput(problemData.t+problemData.dt) * problemData.xiDisc - problemData.cDisc(:,:,1);
[zbDisc, problemData.zbV] = applyVertexAveraging(problemData.g, zbDisc, problemData.averagingOperator, problemData.sysMinValueCorrection);
problemData.changeL2 = norm(problemData.zbDisc - zbDisc, 2);
problemData.zbDisc = zbDisc;

problemData.cDisc(:,:,1) = problemData.rampInput(problemData.t+problemData.dt) * problemData.xiDisc(:,:,1) - problemData.zbDisc;
% Ensure water height doesn't fall below threshold
problemData.cDisc(:,:,1) = correctMinValueExceedanceDisc(problemData.cDisc(:,:,1), problemData.sysMinValueCorrection, nStep, problemData.minTol, problemData.elevTol, problemData.g);
problemData.zbDisc = problemData.rampInput(problemData.t+problemData.dt) * problemData.xiDisc - problemData.cDisc(:,:,1); % TODO consistency in case of correction!

problemData.isSubSteppingFinished = nSubStep >= length(problemData.tLvls);
end % function
