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
  matM = problemData.globS - problemData.stab * problemData.globRmu;
else
  dtA = problemData.dt * problemData.A;
  
  cDiscRkRHS = zeros(K * N, 1);
  for i = 1 : nSubStep - 1
    cDiscRkRHS = cDiscRkRHS + dtA(nSubStep, i) .* problemData.cDiscRK{i};
  end

  matLbar = - problemData.globG{1} - problemData.globG{2} + problemData.stab * problemData.globRphi;
  matL = problemData.globMphi + dtA(nSubStep, nSubStep) .* matLbar;
  vecBphi = problemData.globH - problemData.globFphiIn;
  vecQ =  problemData.globMcDisc + dtA(nSubStep, nSubStep) * vecBphi + cDiscRkRHS;
  matMbar = problemData.globS - problemData.stab * problemData.globRmu;
  matM = dtA(nSubStep, nSubStep) .* matMbar;
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
if problemData.isBlockSolve
  matLinv = blkinv(matL, problemData.blockSolveSize * N);
  LinvQ = matLinv * vecQ;
  LinvM = matLinv * matM;
else
  LinvQ = matL \ vecQ;
  LinvM = matL \ matM;
end % if

%% Solving global system for lambda
matN = - problemData.stab * problemData.globT - problemData.globKmuOut ;
matP = problemData.globP;

lambdaDisc = (-matN * LinvM + matP) \ (problemData.globKmuIn - matN * LinvQ);

%% Reconstructing local solutions from updated lambda
problemData.cDisc = LinvQ - LinvM * lambdaDisc;

if ~problemData.isStationary
  problemData.cDiscRK{nSubStep} = vecBphi - matLbar * problemData.cDisc - matMbar * lambdaDisc;
end % if

problemData.cDisc = reshape(problemData.cDisc, N, K)';
end % function
