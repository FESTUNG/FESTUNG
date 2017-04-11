% Preprocessing of the Runge-Kutta step.

%===============================================================================
%> @file sweVert/preprocessSubStep.m
%>
%> @brief Preprocessing of the Runge-Kutta step.
%===============================================================================
%>
%> @brief Preprocessing of the Runge-Kutta step.
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
%> This routine is executed first in each loop iteration.
%> It takes care of the assembly of time-dependent matrices and right hand
%> side vectors. Furthermore, it evaluates the boundary conditions for the 
%> current Runge-Kutta step.
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures (either filled
%>                      with initial data or the solution from the previous
%>                      loop iteration), as provided by configureProblem()  
%>                      and preprocessProblem(). @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%> @param  nSubStep     The current iteration number of the substepping.
%>
%> @retval problemData  The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
function problemData = preprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
t = problemData.t(nSubStep);

%% L2-projections of algebraic coefficients and right hand side.
DDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) c(t,x1,x2), problemData.N, problemData.qOrd, ...
                       problemData.globM, problemData.basesOnQuad2D), problemData.DCont, 'UniformOutput', false);
problemData.globLu = reshape(projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) problemData.fuCont(t,x1,x2), ...
                             problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D).', [], 1);
problemData.globLh = reshape(projectFuncCont2DataDisc1D(problemData.g.g1D, @(x1) problemData.fhCont(t,x1), ...
                             problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D).', [], 1);

%% Compute depth averaged velocity
barU1Disc = { zeros(problemData.g.g1D.numT, problemData.barN), zeros(problemData.g.g1D.numT, problemData.barN) };
for s = 1 : 2
  for j = 1 : problemData.barN
    i = mapTensorProductIndex(j, 1);
    barU1Disc{s}(:, j) = problemData.g.g1D.markT2DT.' * (problemData.cDiscRK{nSubStep, 2}(:,i) .* problemData.g.J0T{s}(:,2,2));
  end % for j
end % for s

%% Assembly of time-dependent global matrices.
% Advection and diffusion element integrals in momentum equation
problemData.globEP = assembleMatElemTetraDphiPhiFuncDisc(problemData.g, problemData.hatG, problemData.cDiscRK{nSubStep, 2});
problemData.globGR = assembleMatElemTetraDphiPhiFuncDisc(problemData.g, problemData.hatG, DDisc);

% Advection and diffusion interior edge integrals in momentum equation
globR = assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);
problemData.globGR = cellfun(@(g,r) g - r, problemData.globGR, globR, 'UniformOutput', false);
problemData.globEP = cellfun(@(e,p) e - p, problemData.globEP, assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint | (problemData.g.markE0TbdrU & problemData.g.markE0TbdrRiem), problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDiscRK{nSubStep, 2}), 'UniformOutput', false);

% Interior edge flux in continuity and free surface equation
problemData.tildeGlobP = problemData.fn_assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, problemData.g.markE0Tint | (problemData.g.markE0TbdrH & problemData.g.markE0TbdrRiem), problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);
problemData.barGlobP = problemData.fn_assembleMatV0T1DPhiPhiFuncDiscNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0Tint | (problemData.g.g1D.markV0TbdrH & problemData.g.g1D.markV0TbdrRiem), problemData.barHatPdiag, problemData.barHatPoffdiag);

% Advection element integral in free surface equation
problemData.barGlobG = problemData.fn_assembleMatElem1DDphiPhiFuncDiscHeight(barU1Disc, problemData.hSmoothQ0T1D, problemData.barHatG);

%% Assembly of boundary contributions.
hDCont = @(x1) problemData.hDCont(t,x1);
u1DCont = @(x1,x2) problemData.u1DCont(t,x1,x2);
u2DCont = @(x1,x2) problemData.u2DCont(t,x1,x2);
qDCont = { @(x1,x2) problemData.q1DCont(t,x1,x2); @(x1,x2) problemData.q2DCont(t,x1,x2) };
uhDCont = @(x1) problemData.uhDCont(t, x1);

% Diffusion boundary terms in momentum equation
globRbdr = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrQ, problemData.hatRdiag, DDisc);
problemData.globGR = cellfun(@(g,r) g - r, problemData.globGR, globRbdr, 'UniformOutput', false);
globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrQ, qDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = globJtmp{1} + globJtmp{2};

% Advection boundary terms in momentum equation
globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrH & ~problemData.g.markE0TbdrRiem, @(x1,x2) hDCont(x1), problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + problemData.gConst * globJtmp{1};
problemData.globEP = cellfun(@(e,p) e - p, problemData.globEP, assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & (~problemData.g.markE0TbdrU | problemData.g.markE0TbdrRiem), problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2}), 'UniformOutput', false);

globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU & ~(problemData.g.markE0TbdrCoupling | problemData.g.markE0TbdrRiem), @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + globJtmp{1};
problemData.globSbdr = assembleMatEdgeTetraPhiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU & ~(problemData.g.markE0TbdrCoupling | problemData.g.markE0TbdrRiem), u1DCont, problemData.qOrd, problemData.hatQPerQuad);

% Boundary terms for horizontal velocity in continuity and flux equation
problemData.globJu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrCoupling, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU & ~(problemData.g.markE0TbdrCoupling | problemData.g.markE0TbdrRiem), u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJCont = globJtmp{1};

% Dirichlet boundary for vertical velocity in continuity equation
globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrW & ~problemData.g.markE0TbdrCoupling, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJCont = problemData.globJCont + globJtmp{2};

% Boundary terms in free surface equation
problemData.barGlobPbdr = problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0Tbdr & (~problemData.g.g1D.markV0TbdrUH | problemData.g.g1D.markV0TbdrRiem), problemData.barHatPdiag);
problemData.barGlobJuh = zeros(problemData.g.g1D.numT * problemData.barN, 1);
for n = 1 : 2
  markV0TbdrUH = problemData.g.g1D.markV0TbdrUH(:, n) & ~problemData.g.g1D.markV0TbdrRiem(:, n);
  markV0TbdrUHrep = logical(kron(markV0TbdrUH, true(problemData.barN, 1)));
  x1V0T = problemData.g.g1D.coordV0T(markV0TbdrUH, n, 1);
  problemData.barGlobJuh(markV0TbdrUHrep) = ( uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrUH, n) ) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n)';
end % for n

%% Assembly of boundary contributions with Riemann solver
% Advection boundary terms in momentum equation
globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrH & problemData.g.markE0TbdrRiem, @(x1,x2) hDCont(x1), problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + 0.5 * problemData.gConst * globJtmp{1};

globJtmp = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU & problemData.g.markE0TbdrRiem & ~problemData.g.markE0TbdrCoupling, @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + 0.5 * globJtmp{1};

% Boundary terms for horizontal velocity in continuity equation
globJtmp = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrU & problemData.g.markE0TbdrRiem, @(x1,x2) u1DCont(x1,x2) .* hDCont(x1), problemData.hSmoothV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJCont = problemData.globJCont + 0.5 * globJtmp;

% Boundary terms in free surface equation
for n = 1 : 2
  markV0TbdrRiem = problemData.g.g1D.markV0TbdrUH(:, n) & problemData.g.g1D.markV0TbdrRiem(:, n);
  markV0TbdrRiemRep = logical(kron(markV0TbdrRiem, true(problemData.barN, 1)));
  x1V0T = problemData.g.g1D.coordV0T(markV0TbdrRiem, n, 1);
  problemData.barGlobJuh(markV0TbdrRiemRep) = 0.5 * ( uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrRiem, n) ) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n)';
end % for

%% Assembly of jump term in Lax-Friedrichs solver
problemData = problemData.fn_assembleJumpTerms(problemData, nSubStep);
end % function
