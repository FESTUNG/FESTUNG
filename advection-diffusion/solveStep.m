% Second step of the four-part algorithm in the main loop. Computes the
% solution at the next time step.

%===============================================================================
%> @file advection-diffusion/solveStep.m
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>        Computes the solution at the next time step.
%===============================================================================
%>
%> @brief Second step of the four-part algorithm in the main loop.
%>        Computes the solution at the next time step.
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
%> This routine is executed second in each loop iteration.
%> It assembles the global system and solves it using the backslash-operator.
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
function problemData = solveStep(problemData, nStep) %#ok<INUSD>
K = problemData.K;
N = problemData.N;
%% Building and solving the system.

if problemData.reduceAnsatz
    smallNumberAnsatz = nchoosek(problemData.p + 1, problemData.p - 1);
else % if reduceAnsatz
    smallNumberAnsatz = nchoosek(problemData.p + 2, problemData.p);
end % if reduceAnsatz


if problemData.implicitEuler
    
    sysA = [                                                problemData.globM,                                                  sparse(K*N,K*N), -problemData.globH{1}+problemData.globQ{1}+problemData.globQN{1};
        sparse(K*N,K*N),                                                problemData.globM, -problemData.globH{2}+problemData.globQ{2}+problemData.globQN{2};
        -problemData.globG{1}+problemData.globR{1}+problemData.globRD{1}, -problemData.globG{2}+problemData.globR{2}+problemData.globRD{2},                             problemData.globS+problemData.globSD  - problemData.globGadv{1} - problemData.globGadv{2} + problemData.globRadv];
    sysV = [-problemData.globJD{1}; -problemData.globJD{2}; problemData.globKD-problemData.globKN+problemData.globL];
        
    markToCalculate = [ logical(kron(ones(problemData.g.numT, 1), [ones(smallNumberAnsatz, 1); zeros(problemData.N - smallNumberAnsatz, 1)])) ;
        logical(kron(ones(problemData.g.numT, 1), [ones(smallNumberAnsatz, 1); zeros(problemData.N - smallNumberAnsatz, 1)])) ;
        logical(kron(ones(problemData.g.numT, 1), [ones(smallNumberAnsatz, 1); ones(problemData.N - smallNumberAnsatz, 1)])) ] ;
    
    problemData.sysY(markToCalculate) = (problemData.sysW(markToCalculate,markToCalculate) + problemData.tau*sysA(markToCalculate,markToCalculate)) \ ...
        (problemData.sysW(markToCalculate,markToCalculate)*problemData.sysY(markToCalculate) + problemData.tau*sysV(markToCalculate));

else % if implicitEuler
    
    markToCalculate = logical(kron(ones(problemData.g.numT, 1), [ones(smallNumberAnsatz, 1); zeros(problemData.N - smallNumberAnsatz, 1)]));
    
    sysU = problemData.sysY(2*K*N+1 : 3*K*N);
    sysQ1 = zeros(size(sysU));
    sysQ2 = zeros(size(sysU));
    sysQ1(markToCalculate) = problemData.globM(markToCalculate,markToCalculate) \ ( (problemData.globH{1}(markToCalculate,:)-problemData.globQ{1}(markToCalculate,:)-problemData.globQN{1}(markToCalculate,:)) * sysU - problemData.globJD{1}(markToCalculate) );
    sysQ2(markToCalculate) = problemData.globM(markToCalculate,markToCalculate) \ ( (problemData.globH{2}(markToCalculate,:)-problemData.globQ{2}(markToCalculate,:)-problemData.globQN{2}(markToCalculate,:)) * sysU - problemData.globJD{2}(markToCalculate) );
    
    sysUdot = problemData.globM \ ...
        ( problemData.globKD-problemData.globKN+problemData.globL ...
        + ( problemData.globG{1}(:,markToCalculate)-problemData.globR{1}(:,markToCalculate)-problemData.globRD{1}(:,markToCalculate) ) * sysQ1(markToCalculate) ...
        + ( problemData.globG{2}(:,markToCalculate)-problemData.globR{2}(:,markToCalculate)-problemData.globRD{2}(:,markToCalculate) ) * sysQ2(markToCalculate) ...
        - (problemData.globS+problemData.globSD  - problemData.globGadv{1} - problemData.globGadv{2} + problemData.globRadv) * sysU );
    
    % Apply slope limiting to time derivative
    if problemData.isSlopeLim
        cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(sysUdot, [N K])', problemData.globM, problemData.globMDiscTaylor);
        cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim);
        cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
        sysUdot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
    end % if
    
    sysU = sysU + problemData.tau * sysUdot;
    
    % Limiting the solution
    if problemData.isSlopeLim
        % Evaluate boundary condition at new time level
        cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont(problemData.t, x1, x2));
        sysU = reshape(applySlopeLimiterDisc(problemData.g, reshape(sysU, [N K])', problemData.g.markV0TbdrD, ...
            cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim)', [K*N 1]);
    end % if
    
    problemData.sysY = [sysQ1; sysQ2; sysU];

end % if implicit Euler

end