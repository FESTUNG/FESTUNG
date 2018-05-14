% Compute the solution of the current Runge-Kutta stage.

%===============================================================================
%> @file
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
%> This routine is executed second in each loop iteration.
%> It assembles the global system, computes the discrete time derivative
%> and applies slope limiting to it (i.e., applies "selective mass lumping"
%> as described in @ref RAWFK2016). This is used to compute the solution at
%> the next Runge-Kutta level, which then is slope-limited itself.
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

% Building the system
sysA = -problemData.globG{1} - problemData.globG{2} + problemData.globR;
sysV = problemData.globL - problemData.globKD - problemData.globKN;

% Computing the discrete time derivative
cDiscDot = problemData.globM \ (sysV - sysA * problemData.cDiscRK{nSubStep});

% Apply slope limiting to time derivative
if problemData.isSlopeLim
  cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot, [N K])', problemData.globM, problemData.globMDiscTaylor);
  cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
  cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
  cDiscDot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
end % if

% Compute next step
problemData.cDiscRK{nSubStep + 1} = problemData.omega(nSubStep) * problemData.cDiscRK{1} + (1 - problemData.omega(nSubStep)) * (problemData.cDiscRK{nSubStep} + problemData.tau * cDiscDot);

% Limiting the solution
if problemData.isSlopeLim
  % Evaluate boundary condition at new time level
  tBC = getdefault(problemData.t, nSubStep + 1, problemData.t(1) + problemData.tau);
  cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont(tBC, x1, x2));
  problemData.cDiscRK{nSubStep + 1} = reshape(applySlopeLimiterDisc(problemData.g, reshape(problemData.cDiscRK{nSubStep + 1}, [N K])', problemData.g.markV0TbdrD, ...
                                      cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim)', [K*N 1]);
end % if
end % function
