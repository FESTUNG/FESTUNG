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
%> parameter <tt>problemData.isSubSteppingFinished</tt> becomes 
%> <tt>true</tt>.
%> These three steps are:
%>
%>  1. @link swe_2dv/preprocessSubStep.m @endlink
%>  2. @link swe_2dv/solveSubStep.m @endlink
%>  3. @link swe_2dv/postprocessSubStep.m @endlink
%> 
%> This routine is executed second in each loop iteration.
%> It first solves for @f$\vec{q}_h^{(n)}@f$ and @f$(u^2)^{(n)}@f$ and uses those
%> to update @f$(u^1)^{(n+1)}@f$ and @f$h^{(n+1)}@f$, with @f$n@f$ being the previous
%> time level.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      @link swe_2dv/configureProblem.m @endlink and 
%>                      @link swe_2dv/preprocessProblem.m @endlink. 
%%>                     @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with the new solution
%>                      for this Runge-Kutta stage. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2018 Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
%> @author Balthasar Reuter, 2018
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

if problemData.isReducedOrder
  halfOrderN = (floor(problemData.p / 2) + 1).^2;
  markHalfOrder = logical(kron(ones(problemData.g.numT, 1), [ones(halfOrderN, 1); zeros(problemData.N - halfOrderN, 1)]));
  markSetToZero = logical(kron(ones(problemData.g.numT, 1), [zeros(halfOrderN, 1); ones(problemData.N - halfOrderN, 1)]));
end % if

%% Solve for next time level
% Flux variables
qSys = cell(2,1);
for m = 1 : 2
  if problemData.isReducedOrder
    qSys{m} = zeros(problemData.g.numT * problemData.N,1);
    qSys{m}(markHalfOrder) = problemData.globM(markHalfOrder,markHalfOrder) \ (...
        -problemData.globJuFlux{m}(markHalfOrder) - problemData.globJuCoupling{m}(markHalfOrder) + ...
        problemData.globHQ{m}(markHalfOrder,markHalfOrder) * cSys{2}(markHalfOrder));
  else
    qSys{m} = problemData.globM \ (-problemData.globJuFlux{m} - problemData.globJuCoupling{m} + ...
                                    problemData.globHQ{m} * cSys{2});
  end % if
end % for m

% Vertical velocity
cSys{3} = problemData.globHQup \ (problemData.globJu{1} + problemData.globJu{2} + problemData.globVeeJu + problemData.globVeeJh + ...
                                   0.5 * (problemData.globVeeJuRiem + problemData.globVeeJhRiem + problemData.globJuhRiem) + ...
                                   problemData.globJuCoupling{1} + problemData.globJwCoupling + problemData.globKh + ...
                                  (problemData.globHQavg + problemData.globVeeP + problemData.globVeePbdr) * cSys{2});
                                
% Horizontal velocity
if problemData.isReducedOrder
  rhs = problemData.globLu - problemData.globLzBot + problemData.globM \ ( ...
                -problemData.globJ - problemData.globJuuCoupling{1} - problemData.globJuuCoupling{2} - problemData.globKu + ...
                problemData.globVeeHQ * cSys{1} + problemData.globEP{1} * cSys{2} + problemData.globEP{2} * cSys{3} + ...
                problemData.globGR{1} * qSys{1} + problemData.globGR{2} * qSys{2} );
  rhs(markSetToZero) = 0;

  cSys{2} = problemData.omega(nSubStep) * reshape(problemData.cDiscRK{1, 2}.', [], 1) + ...
            (1 - problemData.omega(nSubStep)) * ( cSys{2} + problemData.tau * rhs );
else
  cSys{2} = problemData.omega(nSubStep) * reshape(problemData.cDiscRK{1, 2}.', [], 1) + ...
           (1 - problemData.omega(nSubStep)) * ( cSys{2} + problemData.tau * ( ...
             problemData.globLu - problemData.globLzBot + problemData.globM \ ( ...
               -problemData.globJ - problemData.globJuuCoupling{1} - problemData.globJuuCoupling{2} - problemData.globKu + ...
               problemData.globVeeHQ * cSys{1} + problemData.globEP{1} * cSys{2} + problemData.globEP{2} * cSys{3} + ...
               problemData.globGR{1} * qSys{1} + problemData.globGR{2} * qSys{2} ) ) );
end % if

% Water height
cSys{1} = problemData.omega(nSubStep) * reshape(problemData.cDiscRK{1, 1}.', [], 1) + ...
          (1 - problemData.omega(nSubStep)) * ( cSys{1} + problemData.tau * ( ...
            problemData.globLh + problemData.globBarM \ ( -problemData.globBarJuh - 0.5 * problemData.globBarJuhRiem - ...
              problemData.globBarJu - 0.5 * problemData.globBarJuRiem - problemData.globBarJh - 0.5 * problemData.globBarJhRiem - ...
              problemData.globBarKh + (problemData.globBarG - problemData.globBarP - problemData.globBarPbdr) * cSys{1}) ) );

%% Convert representation vector to representation matrix
problemData.cDiscRK(nSubStep + 1, :) = cellfun(@(c, cRK) reshape(c, size(cRK, 2), size(cRK, 1)).', ...
                                        cSys, problemData.cDiscRK(nSubStep, :), 'UniformOutput', false);
end % function
