% Compute the solution of the current Runge-Kutta stage.

%===============================================================================
%> @file advection/solveSubStep.m
%>
%> @brief Compute the solution of the current Runge-Kutta stage.
%===============================================================================
%>
%> @brief Compute the solution of the current Runge-Kutta stage.
%>
%> The routine iterateSubSteps() repeatedly executes three steps until the
%> parameter <code>problemData.isSubSteppingFinished</code> becomes
%> <code>true</code>.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%>
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop.
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with the new solution
%>                      for this Runge-Kutta stage. @f$[\text{struct}]@f$
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
function problemData = solveSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
K = problemData.K;
Kedge = problemData.g.numE;
N = problemData.N;
stab = problemData.stab;

%% Actual HDG
matL = problemData.globM ./ problemData.dt - problemData.globG{1} - problemData.globG{2} ...
    + stab * problemData.globRphi; % Here goes the time discretization
vecF = problemData.globMcDisc ./ problemData.dt - stab * problemData.globFgamma - stab * problemData.globCd ; % Add here source terms if needed
matM = problemData.globS{1} + problemData.globS{2} - stab * problemData.globRlambda + problemData.globSN{1} + problemData.globSN{2}  + stab * problemData.globRD;

%% Computing local solves
localSolves = mldivide(matL, [vecF matM]);
LinvF = localSolves(:, 1);
LinvM = localSolves(:, 2:end);

%% SolVing global system for lambda
matN = - stab * problemData.globU - problemData.globRgamma ;
matP = problemData.globP;

vecKD = problemData.globKDlambda;

sysMatA = -matN * LinvM + matP;
sysRhs = vecKD - matN * LinvF;

problemData.cDiscLambda = mldivide( sysMatA, sysRhs );

%% Reconstructing local solutions from updated lambda
problemData.cDisc = LinvF - LinvM * problemData.cDiscLambda;
problemData.cDisc = reshape( problemData.cDisc, problemData.N, problemData.g.numT )';

problemData.cDiscLambda = reshape( problemData.cDiscLambda, problemData.Nlambda, problemData.g.numE )';

end % function
