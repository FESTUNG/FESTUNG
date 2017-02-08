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
N = problemData.N;

% %% Solve local systems
% % System construction
% matL = 1.; % Here goes the time discretization
% vecF = 1.; % Here goes the time discretization
% matM = 1.;
% LinvF = 0;
% LinvM = 0;
% 
% %Right hand side for impl. Euler step
% sysLocalB = vecF + (1/problemData.dt) .* problemData.globalM * problemData.cDisc;
% 
% %Possible further modifications for Runge-Kutta time-stepping
% % TODO
% 
% % System solve
% % TODO May I use this syntax? Does Matlab now, which part of the solution
% % shall go into which vector/matrix?
% [LinvF, LinvM] = mldivide(matL, [sysLocalB matM]);
% 
% %% Solve global systems
% matN = 1.;
% matP = 1.;
% 
% matKD = 1.;
% 
% sysA = ( -1. .* matN * LinvM + P );
% sysB = (matKD - matN * LinvF ) ;
% 
% % System solve
% lDiscDot = mldivide(sysA, sysB);
% 
% %Update solution on elements
% cDiscDot = LinvF - LinvM * lDiscDot;
% problemData.cDiscRK{nSubStep + 1} = cDiscDot;
% Testing HDG
matL = problemData.globM ./ problemData.dt - problemData.globG{1} - problemData.globG{2} + problemData.globRphi   ; % Here goes the time discretization
vecF = problemData.globcDiscTime ./ problemData.dt - problemData.globVecFluxDir - problemData.globVecValDir ; % Add here source terms if needed
matM = problemData.globS{1} + problemData.globS{2} - problemData.globRlambda;

localSolves = mldivide(matL, [vecF matM]);
LinvF = localSolves(:, 1);
LinvM = localSolves(:, 2:end);

matN = - problemData.globU - problemData.globRgamma ;
matP = problemData.globP;

vecKD = problemData.globKDlambda;

sysMatA = -matN * LinvM + matP;
sysRhs = vecKD - matN * LinvF;

problemData.cDiscLambda = mldivide( sysMatA, sysRhs );

problemData.cDisc = LinvF - LinvM * problemData.cDiscLambda;
problemData.cDisc = reshape( problemData.cDisc, problemData.N, problemData.g.numT )';

% LinvF = 0;
% LinvM = 0;

%% Standard DG
% 
% % K = problemData.K;
% % N = problemData.N;
% % Building the system
% sysA = -problemData.globG{1} - problemData.globG{2} + problemData.globR;
% sysV = problemData.globL - problemData.globKD - problemData.globKN;
% 
% % Computing the discrete time derivative
% cDiscDot = problemData.globM \ (sysV - sysA * problemData.cDiscRK{nSubStep});
% 
% % Compute next step
% problemData.cDiscRK{nSubStep + 1} = problemData.omega(nSubStep) * problemData.cDiscRK{1} + (1 - problemData.omega(nSubStep)) * (problemData.cDiscRK{nSubStep} + problemData.tau * cDiscDot);

end % function
