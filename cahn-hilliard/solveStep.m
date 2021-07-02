% Second step of the four-part algorithm in the main loop.

%===============================================================================
%> @file ./cahn-hilliard/configureProblem.m
%>
%> @brief Second step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%>
%> This routine is executed second in each loop iteration and is intended to
%> produce the solution at the next step, e.g., at a new time-level.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop.
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next step. @f$[\text{struct}]@f$
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
function problemData = solveStep(problemData, nStep)
% Second step in each loop iteration. Should compute the next step of the
% solution.

% The solveStep can consist of substeps (e.g., in a Runge-Kutta method),
% for which iterateSubSteps() can be used. It is highly recommended to
% store the substep handles in preprocessProblem to reduce call overhead
% drastically.
% problemData.isSubSteppingFinished = false;

fprintf([repmat('#',1,80),'\n','Starting solve step at t = %.3e and tau = %.3e.\n',...
  repmat('#',1,80),'\n'],problemData.t, problemData.tau);

K = problemData.K; N = problemData.N;
newton_version = problemData.newton_version;

% Parameters needed for Newton's method
len = length(problemData.sysY);
time_step_size_tau = problemData.tau;
grid_size_h = problemData.hmax;
initial_guess_Y0 = problemData.sysY;  
matrix_L = [ problemData.epsilon^2 * problemData.globA , - problemData.globM ; ...
             problemData.globM                         , sparse(K*N,K*N)     ];
vector_C = [ problemData.globEold - problemData.globKNCon ; ...
  - problemData.globM * problemData.sysY(1:K*N) - problemData.tau * (problemData.globRHS + problemData.globKNPot)];
tol_abs_L2h = 1e-4;
tol_rel_L2h = 1e-4;
max_num_iterations = 10;
verbose_level = 1;
assembleAmob = @(Y) [];
slopeLimiter = @(Y) slopeLimiterFct(problemData, Y);

switch  newton_version
  case 1
    globAmob = calculateGlobAnonLinear(problemData, problemData.sysY);
    map_Y_to_N_of_Y = @(Y,globAmobNew,beta) assembleNonLinearPart(problemData, Y, globAmob);
    map_Y_to_DN_of_Y = @(Y,globAmobNew,beta) assembleJacobianOfNonLinearPart(problemData, Y, globAmob);
  case 2
    assembleAmob = @(Y) calculateGlobAnonLinear(problemData, Y);
    map_Y_to_N_of_Y = @(Y,globAmob,beta) assembleNonLinearPart(problemData, Y, globAmob);
    map_Y_to_DN_of_Y = @(Y,globAmob,beta) assembleJacobianOfNonLinearPart(problemData, Y, globAmob);
  case 3
    map_Y_to_N_of_Y = @(Y,globAmob,beta) assembleNonLinearPart(problemData, Y);
    map_Y_to_DN_of_Y = @(Y,globAmob,beta) assembleJacobianOfNonLinearPart(problemData, Y);
  case 4
    assembleAmob = @(Y) calculateGlobAnonLinear(problemData, Y);
    map_Y_to_N_of_Y = @(Y,globAmob,beta) beta * assembleNonLinearPart(problemData, Y, globAmob) ...
      + (1-beta) * assembleNonLinearPart(problemData, Y);
    map_Y_to_DN_of_Y = @(Y,globAmob,beta) beta * assembleJacobianOfNonLinearPart(problemData, Y, globAmob) ...
      + (1-beta) * assembleJacobianOfNonLinearPart(problemData, Y);
  otherwise
    error('Newton Version not defined!');
end % switch newton_version

switch problemData.newton_precond
  case 1,  precond = @(Y,globAmob,beta) eye(K*N*2);
  case 2,  precond = @(Y,globAmob,beta) matrix_L;
  case 3,  precond = @(Y,globAmob,beta) matrix_L + map_Y_to_DN_of_Y(Y,globAmob,beta);
  otherwise,  error('Version of Preconditioner not defined');
end % switch problemData.newton_precond

if problemData.zalesak
  yOld = problemData.sysY;
end % if zalesak


if problemData.implicit
  
  if problemData.adaptiveStepping && ~problemData.zalesak
    yOld = problemData.sysY;
  end % if adaptiveStepping
  
  % Solve Cahn-Hillard using Newton's method
  [problemData.sysY, numIter, newton_flag] = ...
     newton_krylov(len, time_step_size_tau, grid_size_h, initial_guess_Y0, matrix_L, ...
       vector_C, map_Y_to_N_of_Y, tol_abs_L2h, tol_rel_L2h, max_num_iterations, verbose_level, ...
       slopeLimiter, precond, newton_version, assembleAmob);

%   if ~problemData.adaptiveStepping
%     assert(newton_flag > -1,'Newton did not converge')
%   end % if ~ adaptiveStepping
  assert(newton_flag ~= 2,'Stationary state reached')
  
  % Stop calculation if stationary state reached for max (adaptive) time step
  if newton_flag == 2 && problemData.tau * problemData.tauIncr > problemData.tauMax
    problemData.isFinished = true;
    fprintf('Stationary state reached at t = %.3e.\n', problemData.t + problemData.tau);
  end

  % Adaptive time stepping
  if problemData.adaptiveStepping && newton_flag < 0
    problemData.tau = problemData.tau * problemData.tauDecr;
    if problemData.tau < problemData.tauMin
      cDisc = reshape(problemData.sysY(1 : K*N), N, K)';
      cLagr = projectDataDisc2DataLagr(cDisc);
      visualizeDataLagr(problemData.g, cLagr, 'c_h', problemData.outputBasename, nStep, problemData.outputTypes);
    end
    assert(problemData.tau >= problemData.tauMin,'Minimum time step size did not bring convergence');
    problemData.sysY = yOld;
    problemData.stepSuccess = false;
  elseif problemData.adaptiveStepping && numIter < max_num_iterations * problemData.newtonStepRatio ...
          && problemData.tau * problemData.tauIncr <= problemData.tauMax
    problemData.t = problemData.t + problemData.tau;
    problemData.tau = problemData.tau * problemData.tauIncr;
    problemData.stepSuccess = true;
  elseif problemData.adaptiveStepping
    problemData.t = problemData.t + problemData.tau;
    problemData.stepSuccess = true;
  end % if adaptiveStepping && convergence/divergence cases

else % explicit time stepping
  
  cVec = problemData.sysY(1:K*N);
  cDisc = reshape(cVec, N, K)';
  
  E = assembleVecElemPhiFuncContPhi(problemData.g, problemData.basesOnQuad, cDisc, problemData.DerivativeOfPsi);
  
  phiVec = problemData.globM \ ( E + problemData.epsilon^2 * problemData.globA * cVec );
  problemData.sysY(K*N+1:end) = phiVec;
  
  globAmob = calculateGlobAnonLinear(problemData, problemData.sysY);
  problemData.sysY(1:K*N) = cVec - problemData.globM \ (problemData.tau * globAmob * phiVec);
  
  assert(norm(problemData.sysY(1:K*N) - cVec) > 1e-8,'Stationary state reached')
  
end % problemData.implicit

if problemData.zalesak
  cDisc = reshape( problemData.sysY(1:K*N), N, K ).';
  
  globFe = assembleCellEdgeTermsLimiting(problemData.g, problemData.basesOnQuad, cDisc, ...
      problemData.mobilityCont, problemData.hatTmobdiag, problemData.hatTmoboffdiag, problemData.hatSdiag, ...
      problemData.hatSoffdiag, problemData.eta, problemData.sigma);
  suppressedFluxes = zeros(K,3);
  phiVec = problemData.sysY(K*N+1:end);
  for n = 1 : 3
    suppressedFluxes(:,n) = (globFe{n} * phiVec) ./ sqrt(2);
  end
  umin = -1;
  umax = +1;
  lowOrderMeans = sqrt(2) * yOld(1:N:K*N);
  [problemData.sysY(1:K*N), flag] = limitFluxFractStepZalesak(problemData.g, umin, umax,...
      problemData.sysY(1:K*N), lowOrderMeans, suppressedFluxes, N, problemData.tau, ...
      1e-7, 1e-7, problemData.zalesakSteps);
   % Store mean values to restore them after slope limiting below.
   sysYmeansAfterFSL = problemData.sysY(1:N:K*N);
end % if zalesak

% Adaptive time stepping
if problemData.adaptiveStepping && (flag == 0)
  problemData.tau = problemData.tau * problemData.tauDecr;
  if problemData.tau < problemData.tauMin
    cDisc = reshape(problemData.sysY(1 : K*N), N, K)';
    cLagr = projectDataDisc2DataLagr(cDisc);
    visualizeDataLagr(problemData.g, cLagr, 'c_h', problemData.outputBasename, nStep, problemData.outputTypes);
  end
  assert(problemData.tau >= problemData.tauMin,'Minimum time step size did not bring convergence');
  problemData.sysY = yOld;
  problemData.stepSuccess = false;
elseif problemData.adaptiveStepping && numIter < max_num_iterations * problemData.newtonStepRatio ...
        && problemData.tau * problemData.tauIncr <= problemData.tauMax
  problemData.t = problemData.t + problemData.tau;
  problemData.tau = problemData.tau * problemData.tauIncr;
  problemData.stepSuccess = true;
elseif problemData.adaptiveStepping
  problemData.t = problemData.t + problemData.tau;
  problemData.stepSuccess = true;
end % if adaptiveStepping && convergence/divergence cases

if problemData.standardLim
  problemData.sysY = slopeLimiterFct(problemData, problemData.sysY);
end % if slope Limiter
% The slope limiter slighlty perturbes the mean values (in the final digit), which
% accumulates over time steps.  Therefore, the mean values are stored after the
% flux limiting and restored after slope limiting.
problemData.sysY(1:N:K*N) = sysYmeansAfterFSL;
  
end % function solveStep

%% Assemble non-linear part
function ret = assembleNonLinearPart(problemData, Y, varargin)

K = problemData.K; N = problemData.N;

cDisc = reshape(Y(1:K*N), N, K)';
phiVec = Y(K*N+1:end);

if isempty(varargin)
  globAmob = calculateGlobAnonLinear(problemData, Y);
  E = assembleVecElemPhiFuncContPhi(problemData.g, problemData.basesOnQuad, cDisc, problemData.ConvexPsi);
elseif length(varargin) == 1
  globAmob = varargin{1};
  E = assembleVecElemPhiFuncContPhi(problemData.g, problemData.basesOnQuad, cDisc, problemData.ConvexPsi);
elseif length(varargin) == 2
  globAmob = varargin{1};
  E = varargin{2};
end % if varargin

ret = [ E ; problemData.tau * globAmob * phiVec ];

end % assembleNonLinearPart

%% Assemble Jacobian of non-linear part
function ret = assembleJacobianOfNonLinearPart(problemData, Y, varargin)

K = problemData.K; N = problemData.N;
cDisc = reshape(Y(1:K*N), N, K)';

if isempty(varargin)
  globAmob = calculateGlobAnonLinear(problemData, Y);
  DerivativeE = assembleMatElemPhiPhiFuncContPhi(problemData.g, problemData.basesOnQuad, ...
    cDisc, problemData.JacobiConvexPsi, problemData.hatQuadPhiPhi);
elseif length(varargin) == 1
  globAmob = varargin{1};
  DerivativeE = assembleMatElemPhiPhiFuncContPhi(problemData.g, problemData.basesOnQuad, ...
    cDisc, problemData.JacobiConvexPsi, problemData.hatQuadPhiPhi);
elseif length(varargin) == 2
  globAmob = varargin{1};
  DerivativeE = varargin{2};
end % if varargin

ret = [ DerivativeE     , sparse(K*N,K*N)            ; ...
        sparse(K*N,K*N) , problemData.tau * globAmob ];

end % function assembleJacobianOfNonLinearPart

%% Assembly of non-linear diffusion matrix
function ret = calculateGlobAnonLinear(problemData, Y)

K = problemData.K; N = problemData.N;
cDisc = reshape( Y(1:K*N), N, K ).';

if problemData.isFluxLim
  globFe = assembleCellEdgeTermsLimiting(problemData.g, problemData.basesOnQuad, cDisc, ...
      problemData.mobilityCont, problemData.hatTmobdiag, problemData.hatTmoboffdiag, problemData.hatSdiag, ...
      problemData.hatSoffdiag, problemData.eta, problemData.sigma);
%   if problemData.implicit
%     alphaEface = getFluxLimitingParameter(Y, problemData.gamma, problemData.g.markE0TE0T, ...
%       globFe, K, N);
%   else
    alphaEface = getFluxLimitingParameter(Y, problemData.gamma, problemData.g.markE0TE0T, ...
      globFe, K, N, problemData.g.areaT);
%   end
else
  alphaEface = ones(K,3);
end % if problemData.isFluxLim

ret = assembleMatMobility(problemData.g, problemData.basesOnQuad, cDisc, problemData.mobilityCont, ...
  problemData.hatLmob, problemData.hatTmobdiag, problemData.hatTmoboffdiag, problemData.hatSdiag, ...
  problemData.hatSoffdiag, problemData.eta, problemData.sigma, alphaEface);

end % function calculateAlphas

%% Function for standtad slope limiting
function ret = slopeLimiterFct(problemData, Y)

N = problemData.N;
K = problemData.K;
ret = Y;
% Evaluate boundary condition at new time level
% tBC = getdefault(problemData.t, nSubStep + 1, problemData.t(1) + problemData.tau);
cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) zeros(size(x1)));
ret(1:K*N) = reshape(applySlopeLimiterDisc(problemData.g, reshape(Y(1:K*N), [N K])', problemData.g.markV0TbdrD, ...
                                      cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim)', [K*N 1]);

% if problemData.standardLim
%   N = problemData.N; K = problemData.K;
%   ret(Y(1:N:K*N) < -1/sqrt(2)) = -1/sqrt(2) + 100 * eps;
%   ret(Y(1:N:K*N) > +1/sqrt(2)) = +1/sqrt(2) - 100 * eps;
% end % if initialLimiting

end % functionSlopeLimiter