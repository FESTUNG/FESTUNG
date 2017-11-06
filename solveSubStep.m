% Compute the solution of the current Runge-Kutta stage.

%===============================================================================
%> @file sweVert/solveSubStep.m
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
%% Convert representation matrix to representation vector
cSys = cellfun(@(c) reshape(c.', [], 1), problemData.cDiscRK(nSubStep, :), 'UniformOutput', false);
%% Solve for next time level
% Flux variables
qSys = cell(2,1);
for m = 1 : 2
  qSys{m} = problemData.globM \ (-problemData.globJuBot{m} - problemData.globJuFlux{m} - problemData.globJuCoupling{m} + ...
                                  problemData.globHQ{m} * cSys{2});
end % for m

% Vertical velocity
cSys{3} = problemData.globHQup \ (problemData.globJu{1} + problemData.globJuBot{1} + 0.5 * problemData.globJuhRiem + problemData.globJw{2} + ...
                                   problemData.globJuCoupling{1} + problemData.globJwCoupling + ...
                                   + (problemData.globKh + problemData.globKhRiem) + ...
                                  (problemData.globHQavg + problemData.tildeGlobP + 0.5 * problemData.tildeGlobPRiem) * cSys{2});
                                
% Horizontal velocity
cSys{2} = problemData.omega(nSubStep) * reshape(problemData.cDiscRK{1, 2}.', [], 1) + ...
          (1 - problemData.omega(nSubStep)) * ( cSys{2} + problemData.tau * ( ...
            problemData.globLu - problemData.globLzBot + problemData.globM \ ( ...
              -problemData.globJ - problemData.globJuuCoupling - problemData.globJuwCoupling - ...
              (problemData.globKu + problemData.globKuRiem) + ...
              problemData.tildeGlobHQ * cSys{1} + problemData.globEP{1} * cSys{2} + problemData.globEP{2} * cSys{3} + ...
              problemData.globGR{1} * qSys{1} + problemData.globGR{2} * qSys{2} ) ) );
%               (problemData.globEP{2} - problemData.globSbdr{2} - 0.5 * problemData.globSriem{2} - problemData.globSCoupling) * cSys{3} + ...
            
% Water height
cSys{1} = problemData.omega(nSubStep) * reshape(problemData.cDiscRK{1, 1}.', [], 1) + ...
          (1 - problemData.omega(nSubStep)) * ( cSys{1} + problemData.tau * ( ...
            problemData.globLh + problemData.barGlobM \ ( -problemData.barGlobJuh - 0.5 * problemData.barGlobJuhRiem - ...
              (problemData.barGlobKh + problemData.barGlobKhRiem) + ...
              (problemData.barGlobG - problemData.barGlobP - problemData.barGlobPbdr - 0.5 * problemData.barGlobPRiem) * cSys{1}) ) );
%% Convert representation vector to representation matrix
problemData.cDiscRK(nSubStep + 1, :) = cellfun(@(c, cRK) reshape(c, size(cRK, 2), size(cRK, 1)).', ...
                                        cSys, problemData.cDiscRK(nSubStep, :), 'UniformOutput', false);
end % function
