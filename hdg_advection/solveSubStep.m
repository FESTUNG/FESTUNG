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

if problemData.isStationary
  matL = - problemData.globG{1} - problemData.globG{2} + problemData.stab * problemData.globRphi;
  vecQ = problemData.globH - problemData.globFphiD;
  matM = problemData.globS + problemData.globSout - problemData.stab * problemData.globRmu;
else
  cDiscRkRHS = zeros(K * N, 1);
  for i = 1 : nSubStep - 1
    cDiscRkRHS = cDiscRkRHS + problemData.A(nSubStep, i) .* problemData.cDiscRK{i};
  end

  matLbar = - problemData.globG{1} - problemData.globG{2} + problemData.stab * problemData.globRphi;
  matL = problemData.globMphi ./ problemData.dt + problemData.A(nSubStep, nSubStep) .* matLbar; % Here goes the time discretization
  vecBphi = problemData.globH - problemData.globFphiD;
  vecQ =  problemData.globMcDisc ./ problemData.dt ...
        + problemData.A(nSubStep, nSubStep) .* vecBphi + cDiscRkRHS; % Add here source terms if needed
  matMbar = problemData.globS ...
                        + problemData.globSout ...
                        - problemData.stab * problemData.globRmu;
  matM = problemData.A(nSubStep, nSubStep) .* matMbar;
end % if

%% Computing local solves
% There are two options.
% 1. Invert the block diagonal matrix L locally, i.e. each block is 
% inverted and then  we construct the inverse matrix L^{-1} from these 
% blocks. This is usuall quickest for large matrices and also saves a lot 
% of memory. It may be efficient to invert more than one block at once.
% 2. We invert the whole mass matrix. This may be very slow and memory
% consuming for large matrices (=many elements). I guess it may be faster
% for matrices of moderate size.
if problemData.isTrueLocalSolve
    blockSize = problemData.localSolveBlockSize;
    matLinvLocal = cell( K/blockSize, 1);
    %Invert every block locally
    for iE=1:K/blockSize
        iEs = (iE-1)*N*blockSize + 1; %index Element start
        iEe =  iE*N*blockSize; %index Element end
        matLinvLocal{iE} =  mldivide(matL(iEs:iEe,iEs:iEe), speye(N*blockSize,N*blockSize) );
    end
    %Construct inverse matrix
    matLinv = blkdiag(  matLinvLocal{:} );
    %Solve L x = [vecQ matM]
    localSolves = matLinv * [vecQ matM];
    LinvQ = localSolves(:, 1);
    LinvM = localSolves(:, 2:end);
else
    localSolves = mldivide(matL, [vecQ matM]);
    LinvQ = localSolves(:, 1);
    LinvM = localSolves(:, 2:end);
end % if
%% Solving global system for lambda
matN = - problemData.stab * problemData.globT - problemData.globKmuOut ;
matP = problemData.globP;

vecR = problemData.globKmuD;

sysMatA = -matN * LinvM + matP;
sysRhs = vecR - matN * LinvQ;

problemData.lambdaDisc = mldivide( sysMatA, sysRhs );

%% Reconstructing local solutions from updated lambda
problemData.cDisc = LinvQ - LinvM * problemData.lambdaDisc;
problemData.cDisc = reshape( problemData.cDisc, problemData.N, problemData.g.numT )';
problemData.lambdaDisc = reshape( problemData.lambdaDisc, problemData.Nmu, problemData.g.numE )';
end % function
