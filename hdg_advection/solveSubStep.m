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

if (problemData.isInTesting == true)
    %% Actual HDG for testing with dt
    matL = problemData.globM ./ problemData.dt - problemData.globG{1} - problemData.globG{2} ...
        + stab * problemData.globRphi; % Here goes the time discretization
    vecF = problemData.globMcDisc ./ problemData.dt - stab * problemData.globFgamma - stab * problemData.globCd ; % Add here source terms if needed
    matM = problemData.globS{1} + problemData.globS{2} - stab * problemData.globRlambda;
    
    %% Computing local solves
    localSolves = mldivide(matL, [vecF matM]);
    LinvF = localSolves(:, 1);
    LinvM = localSolves(:, 2:end);
    
    %% SolVing global system for lambda
    matN = - stab * problemData.globU - problemData.globRgamma ;
    matP = problemData.globP;
    
    vecKD = problemData.globKDlambda;
    
    resC11 = matL * problemData.cDiscReshaped;
    resC12 = matM * problemData.cDiscLambdaReshaped;
    resC = (resC11 + resC12) - vecF;
    
    
    resEdge11 = (- problemData.globG{1} - problemData.globG{2})* problemData.cDiscReshaped;
    resEdge12 = problemData.globS{1} * problemData.cDiscLambdaReshaped;
    resEdge13 = problemData.globS{2} * problemData.cDiscLambdaReshaped;
    
    resEdge21 = stab * problemData.globFgamma;
    
    resEdge = ( resEdge11 + resEdge12  + resEdge13) + resEdge21;
    
    resHyb1 = matN * problemData.cDiscReshaped + matP * problemData.cDiscLambdaReshaped;
    resHyb2 = vecKD;
    resHyb = resHyb1 - resHyb2;
    
    r11 = problemData.globRphi * problemData.cDiscReshaped;
    r12 = problemData.globRlambda * problemData.cDiscLambdaReshaped;
    
%     r21 = problemData.globMcDisc ./ problemData.dt;
%     r22 = problemData.globFgamma;
    r21 = 0;
    r22 = 0;
    r23 = problemData.globCd;
    
    diffr1r2 = r11 - r12 - r21 + r22 + r23;
    
    sysMatA = -matN * LinvM + matP;
    sysRhs = vecKD - matN * LinvF;
    
    cDiscLambdaOld = problemData.cDiscLambda;
    cDiscLambdaOldReshaped = reshape( cDiscLambdaOld', size(problemData.globP, 1), 1 );
    cDiscLambdaUnused = mldivide( sysMatA, sysRhs );
    
    %% Testing compute lambda given exact cDisc
    lambdaFromC = mldivide( matP, vecKD - matN*problemData.cDiscReshaped );
    %% Debugging
    
    % Inserting 'approximated' cDisc and cDiscLambda to compare residuals
    r1d = problemData.globRphi * problemData.cDiscReshaped - problemData.globRlambda * cDiscLambdaUnused;
    r2d = - stab * problemData.globFgamma - problemData.globCd;
    diffr1r2d = r1d - r2d;
    
    %Inserting exact lambda and reconstruct c
    warning('Inserting exact solution for lambda')
    problemData.cDiscLambda = projectFuncCont2FaceDataDisc(problemData.g, @(x1, x2) problemData.getRGSol(problemData.t+problemData.dt, x1, x2), problemData.p*2, problemData.hatMlambda, problemData.basesOnGamma);
    problemData.cDiscLambda = reshape( problemData.cDiscLambda', problemData.Nlambda*problemData.g.numE, 1 );
    CFromLambda = mldivide( matL, vecF - matM*problemData.cDiscLambda );
    
    % Inserting 'exact' cDisc and cDiscLambda to compare residuals
    % r1 = problemData.globRphi * reshape( problemData.cDisc', size(problemData.globM, 1), 1 ) - problemData.globRlambda * problemData.cDiscLambda;
    % r2 = - problemData.globFgamma - problemData.globCd;
    % diffr1r2 = r1 - r2;
   
    %Compute residual
    % problemData.t+problemData.dt
    % matL*reshape( problemData.cDisc', size(problemData.globM, 1), 1 )
    % - matM * problemData.cDiscLambda
    res = vecF- matL*reshape( problemData.cDisc', size(problemData.globM, 1), 1 ) - matM * problemData.cDiscLambda;
    
    %% Reconstructing local solutions from updated lambda
    problemData.cDisc = LinvF - LinvM * problemData.cDiscLambda;
    problemData.cDisc = reshape( problemData.cDisc, problemData.N, problemData.g.numT )';
    
    problemData.cDiscLambda = reshape( problemData.cDiscLambda, problemData.Nlambda, problemData.g.numE )';
else
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
end


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
