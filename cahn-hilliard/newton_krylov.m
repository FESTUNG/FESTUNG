%> Author: Florian Frank
%>
%> THIS CODE SHOULD NOT BE PUBLISHED.
%>
%>   @brief Newton's method.
%>
%>   We search the root of [L*Y + N(Y) + C] = 0, where L*Y is the linear,
%>   N(Y) the nonlinear, and C the constant part of the system matrix of
%>   the current time step (L, DN are matrices, N, C, Y are vectors).
%>   Within a Newton iteration, L and C stay constant
%>   and only N(Y) and its Jacobian DN(Y) have to be updated.  We call the
%>   iterated solutions Y_old and Y_new.
%>   The following storage is required:
%>   - Y_old and Y_new, the first Y_old will be given by initial_guess_Y0.
%>   - residual_old = residual(Y_old) and residual_new = residual(Y_new), which
%>     the negative right-hand sides of the Newton system of equations below.
%>   - initial_norm_residual_L2h, new_norm_residual_L2h, the L2h norm of the residuals.
%>   There is a flag indicating the status of the algorithm.  The algorithm will
%>   will terminate as soon as the flag is set to indicate a success or failure
%>   criterion.
%>
%>   Algorithm:
%>
%>   (1) Initialization:
%>       - Initialize Y_old by initial_guess_Y0.
%>       - Set s = 0, where s is the iteration index (work directly with output argument).
%>       - Compute residual := residual(Y_old) including N(Y_old).
%>       - Compute the norm of residual and store it as initial_norm_residual_L2h
%>         and new_norm_residual_L2h (we need to store the norm of the initial
%>         residual for the relative convergence criterion).
%>       - Check if initial (absolute) residual is below tolerance tol_abs (this
%>         happens at stationary state).  If so, exit with success.
%>
%>   (2) Assembly and solving:
%>       - Build the system of equations including DN(Y_old).
%>       - Solve [L + DN(Y_old)]*Delta = L*Y_old + N(Y_old) + C =: residual(Y_old).
%>
%>   (3) Update and residuals (line_search_Armijo):
%>       - Update Y_new = Y_old - lambda*Delta.    old_norm_residual_L2h = new_norm_residual_L2h;
%>         This is done with the Armijo rule.  If lambda = 1 yields a decrease
%>         in the residual, then lambda = 1 is used.  Otherwise, lambda will be
%>         decreased incrementally until monotonicity is reached.  A monotonicity test
%>         in (4) is not necessary.
%>       - Compute residual = residual(Y_new) (will be right-hand side in next iteration).
%>       - Compute new_norm_residual_L2h = norm(residual).
%>
%>   (4) Tests:
%>       - Convergence: If new_norm_residual_L2h < tolerance_abs or
%>                         new_norm_residual_L2h < tolerance_rel*initial_norm_residual_L2h
%>                      then exit with success.
%>       - Max. iterations: If s == max_iterations then exit with failure.
%>
%>   (5) Schedule and start next iteration.
%>       - Increment s.
%>       - norm_residual_old = norm_residual_new.
%>       - Y_old = Y_new.
%>       - Go to (2).
%>
%> @retval solution_Y               The approximates solution of the nonlinear system (must be allocated but has not to be initialized).
%> @retval num_iterations           Number of iterations the Newton solver took.
%> @retval flag                 Flag: 1:
%> @param  len              The lenth of vector Y.
%> @param  time_step_size_tau       The error tolerances for Newton will be scaled by the time step size.
%> @param  grid_size_h              The mesh size of the underlying regular grid (required for computing the L2h norm).
%> @param  initial_guess_Y0         Initial guess for the Newton iteration (often solution of the previous time step).
%> @param  matrix_L                 Matrix of the linear part of the system.
%> @param  vector_C                 Constant part of the system.
%> @param  map_Y_to_N_of_Y          Nonlinear part of the system as a function of Y.
%> @param  map_Y_to_DN_of_Y         Jacobian of the nonlinear part as function of Y.
%> @param  tol_abs_L2h              Tolerance for the absolute Newton residual error measured in the <b>L2h norm</b>.
%> @param  tol_rel_L2h              Tolerance for the relative Newton residual error measured in the <b>L2h norm</b>.
%> @param  max_num_iterations       Maximum number of allowed Newton iterations; if the solution did not converge when reached, then exit with failure.
%> @param  verbose_level            Silent (0) or print info (1).

function [Y_new, num_iterations_s, flag] = ...
  newton_krylov(len,  time_step_size_tau, grid_size_h, initial_guess_Y0, matrix_L, vector_C, map_Y_to_N_of_Y, ...
  tol_abs_L2h, tol_rel_L2h, max_num_iterations, verbose_level, slopeLimiter, precond, newton_version, assembleAmob)

% Return flags.
SUCCESS_ABSOLUTE_TOLERANCE = 1;
SUCCESS_ABSOLUTE_TOLERANCE_ZERO_STEPS = 2;
SUCCESS_RELATIVE_TOLERANCE = 3;
RUNNING = 0;
FAILURE_MAXIMUM_ITERATIONS = -1;
FAILURE_ZERO_STEP_SIZE_REDUCTION = -2;

% Assertions.
assert(isequal(size(initial_guess_Y0), [len,   1]),   'Dimension mismatch.');
assert(isequal(size(matrix_L),         [len,   len]), 'Dimension mismatch.');
assert(isequal(size(vector_C),         [len,   1]),   'Dimension mismatch.');
assert(tol_abs_L2h        > 0, 'Argument tol_abs_Lh2 must be positive.');
assert(tol_rel_L2h        > 0, 'Argument tol_rel_Lh2 must be positive.');
assert(max_num_iterations > 0, 'Argument max_num_iterations must be positive.');
assert(verbose_level == 0 || verbose_level == 1, 'Verbose level must be 0 or 1.');
assert(isa(slopeLimiter, 'function_handle'), 'slope Limiter must be given as function handle');
assert(isa(precond, 'function_handle'),      'preconditioner must be given as function handle');
assert(newton_version == 1 || newton_version == 2 || newton_version == 3 || newton_version == 4, ...
  'Newton Version has to be an integer between 0 and 5');
assert(isa(assembleAmob, 'function_handle'), 'Assembly of Amob must be given as function handle');

% Define parameters for Newton methods
gamma = 0.9; % Gamma parameter in Eisenstadt Walker algorithm
etamax = 0.1; % Initial relative residual for iterative Krylow solver
max_iters_fdgmres = 200; % maximum number of iterations for gmres
max_restarts_fdgmres = 4; % maximum number of restarts for gmres

beta = 1; %TODO

% Initialization:
num_iterations_s = 0;                                                    % This is an output argument.
Y_old = initial_guess_Y0;                                                % Y^{s-1} = Y0.

% Compute residual vector and norm of residual, store the initial residual.
[vector_R_of_Y, globAmob] = ...
  computeResidual(Y_old, matrix_L, map_Y_to_N_of_Y, vector_C, len, newton_version, ...
  [], beta, assembleAmob);                                         % Residual vector R(Y).

old_norm_residual_L2h = grid_size_h^(2/2) * norm(vector_R_of_Y);         % Using the L2h norm here, ||.||_L2h = h^(2/2)*||.||_2.
initial_norm_residual_L2h = old_norm_residual_L2h;                       % Store the initial norm of residual.
% beta = (num_iterations_s < 20);

if verbose_level == 1
  fprintf('Newton solver starts with initial (absolute) residual: %.3e.\n', initial_norm_residual_L2h);
end % if

% Convergence test (absolute criterion applied on initial residual).
if initial_norm_residual_L2h < time_step_size_tau*tol_abs_L2h
  if verbose_level == 1
    fprintf('Newton converges after 0 steps by SUCCESS_ABSOLUTE_TOLERANCE (residual: %.3e).\n', initial_norm_residual_L2h);
    warning('Warning: a stationary state was reached.  If this is not a physical stationary state then reduce ''nlsolve_tol_abs''.')
  end % if
  
  % Set output arguments.
  Y_new = Y_old;
  flag = SUCCESS_ABSOLUTE_TOLERANCE_ZERO_STEPS;
  return
end % if

% Start Newton loop.
flag = RUNNING;                                                          % The return value.
% P = precond(Y_old,globAmob,beta);
% P = @(x) P\x;
while flag == RUNNING
  num_iterations_s = num_iterations_s + 1;
  if verbose_level == 1
    fprintf('Newton iteration %d.\n', num_iterations_s);
  end % if
  
  P = precond(Y_old,globAmob,beta);
  P = @(x) P\x;
  
  % Solve Newton system [L + DN(Y_old)]*Delta = L*Y_old + N(Y_old) + C.
  %P = precond(Y_old);
  [vector_Delta, error, iter] = fdkrylov(vector_R_of_Y,...
    @(Y) computeResidual(Y, matrix_L, map_Y_to_N_of_Y, vector_C, len, newton_version, ...
         globAmob, beta, assembleAmob),...
    Y_old, [abs(etamax) *  old_norm_residual_L2h,max_iters_fdgmres,max_restarts_fdgmres],2,P);
  vector_Delta = -(P(vector_Delta));
  
  if verbose_level == 1
    fprintf([repmat('+',1,80),'\n'...
      'FD-GMRES error: %.3e, etamax: %.3e, [%d;%d] itertations used from [%d;%d]\n',...
      repmat('+',1,80),'\n'],...
      error, abs(etamax),iter(1),iter(2), max_iters_fdgmres,max_restarts_fdgmres);
  end
  % Line search (Armijo rule).  The norm of the new residual is used to
  % test for convergence, the residual vector is the right-hand side for the
  % next Newton iteration.  A monotonicity test is not necessary, as
  % the Armijo rule guarantees monotonicity.
  [Y_new, vector_R_of_Y, new_norm_residual_L2h, flag] = ...
    linesearchArmijo(Y_old, vector_Delta, matrix_L, map_Y_to_N_of_Y, vector_C, grid_size_h, len, verbose_level, ...
      newton_version, globAmob, beta, assembleAmob);
  
  if verbose_level == 1
    fprintf('Newton absolute residual: %.3e, relative residual: %.3e.\n', ...
      new_norm_residual_L2h, new_norm_residual_L2h/initial_norm_residual_L2h);
  end
  
  %   if isFluxLim || nonConstantMob
  %     enhanceValue = calculateEnhancement(Y_new);
  %     vector_R_of_Y = computeResidual(Y_new, matrix_L, map_Y_to_N_of_Y, vector_C, len, enhanceValue);   % Residual vector R(Y).
  %     new_norm_residual_L2h = grid_size_h^(2/2) * norm(vector_R_of_Y);
  %     if verbose_level == 1
  %       fprintf('Newton absolute residual after new evaluation of enhanceValue: %.3e, relative residual: %.3e.\n', ...
  %         new_norm_residual_L2h, new_norm_residual_L2h/initial_norm_residual_L2h);
  %     end
  %   end
  
  % Standard slope limiting
  %     if new_norm_residual_L2h < time_step_size_tau*tol_abs_L2h
%   Y_new = slopeLimiter(Y_new);
%   [vector_R_of_Y, globAmob] = computeResidual(Y_new, matrix_L, map_Y_to_N_of_Y, vector_C, len, newton_version, [], beta, assembleAmob);
%   new_norm_residual_L2h = grid_size_h^(2/2) * norm(vector_R_of_Y);
%   if verbose_level == 1
%     fprintf('Newton absolute residual after standard limiting procedure: %.3e, relative residual: %.3e.\n', ...
%       new_norm_residual_L2h, new_norm_residual_L2h/initial_norm_residual_L2h);
%   end
  %     end
  
  % Convergence test (absolute criterion).
  if new_norm_residual_L2h < time_step_size_tau*tol_abs_L2h
    if verbose_level == 1
      fprintf('Newton converges after %d steps by absolute criterion.  Exiting with SUCCESS_ABSOLUTE_TOLERANCE.\n', num_iterations_s);
      flag = SUCCESS_ABSOLUTE_TOLERANCE;
      continue
    end % if
  end % if
  
  % Convergence test (relative criterion).
  if new_norm_residual_L2h/initial_norm_residual_L2h < time_step_size_tau*tol_rel_L2h
    if verbose_level == 1                                                % Print info every param.num_out_steps_-th step only.
      fprintf('Newton converges after %d steps by relative criterion.  Exiting with flag SUCCESS_RELATIVE_TOLERANCE.\n', num_iterations_s);
      flag = SUCCESS_RELATIVE_TOLERANCE;
      continue
    end % if
  end % if
  
  % Check if the next iteration would exceed the maximum number of allowed iterations.
  if num_iterations_s == max_num_iterations
    if verbose_level == 1
      fprintf('Maximum number of Newton iterates reached.  Exiting with FAILURE_MAXIMUM_ITERATIONS.\n');
      flag = FAILURE_MAXIMUM_ITERATIONS;
      continue
    end % if
  end % if
  
  
  % Initialize next iterate.
  
  stop_tol_L2h = tol_abs_L2h + tol_rel_L2h * new_norm_residual_L2h;
  ratio_norm_residual_L2h = new_norm_residual_L2h/old_norm_residual_L2h;
  if etamax > 0
    etaold = abs(etamax);
    etanew = gamma * ratio_norm_residual_L2h.^2;
    if gamma * etaold * etaold > .1
      etanew = max(etanew,gamma * etaold * etaold);
    end
    etamax = min([etanew,etamax]);
    etamax = min(etamax,max(0.5 * stop_tol_L2h/new_norm_residual_L2h,etamax));
  end
  Y_old = Y_new;
  if new_norm_residual_L2h / old_norm_residual_L2h > 0.9
    beta = 0;
  end
  old_norm_residual_L2h = new_norm_residual_L2h;
end % while

if verbose_level == 1
  disp('Newton solver ends.');
end % if

end % function



%> @brief Armijo line search/damping to find suitable step length.
%>
%> See C.T.Kelley, "Solving Nonlinear Equations with Newton's Method", Alg. 1.1, lines 6-11.
%>
%> @retval vector_X_new           This output argument will be allocated within this function.
%> @retval vector_R_of_X_new      The residual vector of vector_X_new.
%> @retval new_norm_residual_L2h  Norm of the residual measure in the L2h norm.
%> @param verbose_level           See newton().
function [vector_X_new, vector_R_of_X_new, new_norm_residual_L2h, flag] ...
  = linesearchArmijo(vector_X_old, vector_Delta, matrix_L, map_Y_to_N_of_Y, vector_C, grid_size_h, len, verbose_level,...
      newton_version, globAmob, beta, assembleAmob)

% Asserts.
assert(isequal(size(vector_X_old), [len,   1]),   'Dimension mismatch.');
assert(isequal(size(vector_C),     [len,   1]),   'Dimension mismatch.');
assert(isequal(size(matrix_L),     [len,   len]), 'Dimension mismatch.');

% Constants.
sigma = 0.5;                                                             % A magic number which should be between 0.1 and 0.5.
step_length_lambda = 1.0;                                                % This will be adjusted.
alpha = 1E-4;                                                            % Magic number proposed by C. T. Kelley.
% Decide which enhancement is applied
% Compute old residual vector and norm of old residual.
vector_R_of_X_old = computeResidual(vector_X_old, matrix_L, map_Y_to_N_of_Y, vector_C, len, newton_version, ...
         globAmob, beta, assembleAmob);
old_norm_residual_L2h = grid_size_h^(2/2)*norm(vector_R_of_X_old);       % Using the L2h norm here, ||.||_L2h = h^(2/2)*||.||_2.
flag = 0;

% Compute new residual vector and norm of new residual.
vector_X_new = vector_X_old;                                                % Copy.
vector_X_new = vector_X_new - step_length_lambda * vector_Delta;
vector_R_of_X_new = computeResidual(vector_X_new, matrix_L, map_Y_to_N_of_Y, vector_C, len, newton_version, ...
         globAmob, beta, assembleAmob);
new_norm_residual_L2h = grid_size_h^(2/2) * norm(vector_R_of_X_new);     % Using the L2h norm here, ||.||_L2h = h^(2/2)*||.||_2.

%Monotonicity test.
while new_norm_residual_L2h > (1.0 - alpha*step_length_lambda)*old_norm_residual_L2h
  step_length_lambda = step_length_lambda * sigma;                       % Decreasing the step length.
  vector_X_new = vector_X_old;
  vector_X_new = vector_X_new - step_length_lambda * vector_Delta;
  vector_R_of_X_new = computeResidual(vector_X_new, matrix_L, map_Y_to_N_of_Y, vector_C, len, newton_version, ...
         globAmob, beta, assembleAmob);
  new_norm_residual_L2h = grid_size_h^(2/2) * norm(vector_R_of_X_new);   % Using the L2h norm here, ||.||_L2h = h^(2/2)*||.||_2.
end % while

if step_length_lambda < 0.95 && verbose_level == 1
  fprintf('Necessary step size reduction: %.1f.\n', step_length_lambda);
end % if
% Check if step size reduction is small.
if step_length_lambda < 1e-9
  if verbose_level == 1
    fprintf('Step size reduction smaller than 1e-9.  Exiting with FAILURE_ZERO_STEP_SIZE_REDUCTION.\n');
    flag = -2;
  end % if
end % if
end % function

%> @brief Compute the residual vector R(Y) = L*Y + N(Y) + C at position Y.
%>
%> @param len Size of the solution vector.
function [vector_R_of_Y, varargout] = computeResidual(vector_Y, matrix_L, map_Y_to_N_of_Y, vector_C, len, ...
  newton_version, globAmob, beta, assembleAmob)

assert(isequal(size(vector_Y), [len,   1]), 'Dimension mismatch.');
assert(isequal(size(vector_C), [len,   1]), 'Dimension mismatch.');
assert(isequal(size(matrix_L), [len, len]), 'Dimension mismatch.');

if ( newton_version == 2 && isempty(globAmob) ) || newton_version == 4
  globAmob = assembleAmob(vector_Y);
end % if globAmob

if nargout == 2
  varargout{1} = globAmob;
end % if nargout == 2

vector_R_of_Y = matrix_L*vector_Y + map_Y_to_N_of_Y(vector_Y, globAmob, beta) + vector_C;

end % function