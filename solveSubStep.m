% Second step of the three-part algorithm in the iterateSubSteps loop of each
% step in the main loop.

%===============================================================================
%> @file template/solveSubStep.m
%>
%> @brief Second step of the three-part algorithm in the iterateSubSteps loop of
%>        each step in the main loop.
%===============================================================================
%>
%> @brief Second step of the three-part algorithm in the iterateSubSteps loop of
%>        each step in the main loop.
%>
%> The iterateSubSteps loop repeatedly executes three steps until the number of
%> substep iterations equals the order of the underlying Runge-Kutta method.
%> These three steps are:
%>
%>  1. preprocessSubStep()
%>  2. solveSubStep()
%>  3. postprocessSubStep()
%> 
%> This routine is executed second in each substep loop iteration and is
%> intended to produce the solution at the next substep, e.g., at the time-level
%> of the new Runge-Kutta step.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nSubStep     The current iteration number of the iterateSubSteps 
%>                      loop. 
%>
%> @retval problemData  The input struct enriched with solution data at
%>                      the next substep. @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
%>
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
function problemData = solveSubStep(problemData, ~, nSubStep)
K = problemData.K;
N = problemData.N;
p = problemData.p;

qOrd2D = max(2* p,1); 
[Q1, Q2, W] = quadRule2D(qOrd2D); 
numQuad2D = length(W);

mask = logical(problemData.mask * ones(problemData.numSpecies,1));

% Assembly of time-dependent global matrices
globG = assembleMatElemDphiPhiFuncDiscVec(mask, problemData.g, problemData.hatG, problemData.uHDisc, problemData.vHDisc);
globR = assembleMatEdgePhiPhiValUpwind(mask, problemData.g, problemData.g.markE0Tint, problemData.hatRdiagOnQuad, ...
                                       problemData.hatRoffdiagOnQuad,  problemData.vNormalOnQuadEdge, problemData.g.areaE0TbdrNotN);
% Building the system
sysA = -globG{1} - globG{2} + globR;

for species = 1:problemData.numSpecies
  
  % L2 projections of algebraic coefficients
  fDisc  = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.fCont{species}(problemData.timeLvls(nSubStep),x1,x2), ...
                                    2*p, problemData.hatM, problemData.basesOnQuad);

  % Assembly of Dirichlet boundary contributions
  globKD = assembleVecEdgePhiIntFuncContVal(problemData.g, problemData.g.markE0TbdrD, ...
            @(x1,x2) problemData.cDCont{species}(problemData.timeLvls(nSubStep),x1,x2),  problemData.vNormalOnQuadEdge, N, ...
            problemData.basesOnQuad, problemData.g.areaE0TbdrD);

  % Assembly of Neumann boundary contributions
  gNUpwind = @(x1,x2) (problemData.gNCont{species}(problemData.timeLvls(nSubStep),x1,x2) <= 0) .* problemData.gNCont{species}(problemData.timeLvls(nSubStep),x1,x2);
  globKN = assembleVecEdgePhiIntFuncCont(problemData.g, problemData.g.markE0TbdrN, ...
            gNUpwind, N, problemData.basesOnQuad);

  % Assembly of the source contribution
  globL = problemData.globM * reshape(fDisc', K*N, 1) ...
        + problemData.globT * reshape(problemData.reactions{species}(problemData.timeLvls(nSubStep), problemData.g.mapRef2Phy(1, Q1, Q2), problemData.g.mapRef2Phy(2, Q1, Q2), problemData.cQ0T, problemData.cHQ0T).', numQuad2D*K, 1);
  
  % right hand side
  sysV = globL - globKD - globKN;

  numElem = sum(mask);
  advectiveTerm = zeros(K*N,1);
  advectiveTerm(logical(kron(mask,ones(N,1)))) = sysA * reshape(problemData.concentrationDiscRK{species}(mask,:)', [numElem*N 1]);
  % Computing the discrete time derivative
  cDiscDot = problemData.globM \ (sysV - advectiveTerm);

  % Apply slope limiting to time derivative
  if problemData.isSlopeLim{species}
    cDiscDotTaylor = projectDataDisc2DataTaylor(reshape(cDiscDot, [N K])', problemData.globM, problemData.globMDiscTaylor);
    cDiscDotTaylorLim = applySlopeLimiterTaylor(problemData.g, cDiscDotTaylor, problemData.g.markV0TbdrD, NaN(K,3), problemData.basesOnQuad, problemData.typeSlopeLim{species});
    cDiscDotTaylor = reshape(cDiscDotTaylorLim', [K*N 1]) + problemData.globMCorr * reshape((cDiscDotTaylor - cDiscDotTaylorLim)', [K*N 1]);
    cDiscDot = reshape(projectDataTaylor2DataDisc(reshape(cDiscDotTaylor, [N K])', problemData.globM, problemData.globMDiscTaylor)', [K*N 1]);
  end % if

  % Compute next step
  problemData.cDiscRK{species} = problemData.omega(nSubStep) * problemData.cDiscRK0{species} + (1 - problemData.omega(nSubStep)) * (problemData.cDiscRK{species} + problemData.tau * cDiscDot);

  % Limiting the solution
  if problemData.isSlopeLim{species}
    cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont{species}(problemData.timeLvls(nSubStep), x1, x2));
    [cDiscRK, minMaxV0T] = applySlopeLimiterDisc(problemData.g, reshape(problemData.cDiscRK{species}, [N K])', problemData.g.markV0TbdrD, ...
                                                 cDV0T, problemData.globM, problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim{species});
    problemData.cDiscRK{species} = reshape(cDiscRK', [K*N 1]);
    if problemData.isMask
      problemData.mask(:,species) = (max(minMaxV0T{2},[],2) - min(minMaxV0T{1},[],2)) >= problemData.maskTol;
    end % if
  end % if
  problemData.numOperations = problemData.numOperations + numElem;
end % for

end % function
