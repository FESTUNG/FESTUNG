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
K = problemData.g.numT;
barK = problemData.g.g1D.numT;

%% Algebraic right hand sides and boundary conditions at current time.
fuCont = @(x1,x2) problemData.fuCont(t,x1,x2);
fhCont = @(x1) problemData.fhCont(t,x1);

hDCont = @(x1) problemData.hDCont(t,x1);
u1DCont = @(x1,x2) problemData.u1DCont(t,x1,x2);
u2DCont = @(x1,x2) problemData.u2DCont(t,x1,x2);
qDCont = { @(x1,x2) problemData.q1DCont(t,x1,x2); @(x1,x2) problemData.q2DCont(t,x1,x2) };
uhDCont = @(x1) problemData.uhDCont(t, x1);

%% L2-projections of algebraic coefficients and right hand side.
DDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) c(t,x1,x2), problemData.N, problemData.qOrd, ...
                       problemData.globM, problemData.basesOnQuad2D), problemData.DCont, 'UniformOutput', false);
problemData.globLu = reshape(projectFuncCont2DataDiscTetra(problemData.g, fuCont, problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D).', [], 1);
problemData.globLh = reshape(projectFuncCont2DataDisc1D(problemData.g.g1D, fhCont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D).', [], 1);

%% Compute depth integrated velocity
barU1Disc = { zeros(problemData.g.g1D.numT, problemData.barN), zeros(problemData.g.g1D.numT, problemData.barN) };
for s = 1 : 2
  for j = 1 : problemData.barN
    i = mapTensorProductIndex(j, 1);
    barU1Disc{s}(:, j) = problemData.g.g1D.markT2DT.' * (problemData.cDiscRK{nSubStep, 2}(:,i) .* problemData.g.J0T{s}(:,2,2));
  end % for j
end % for s

%% Assembly of time-dependent matrices and vectors in momentum equation.
% Advection element integral (II)
problemData.globEP = assembleMatElemTetraDphiPhiFuncDisc(problemData.g, problemData.hatG, problemData.cDiscRK{nSubStep, 2});

% Advection interior edge integrals (averaging uu, uw) (VI)
globP = assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDiscRK{nSubStep, 2});
problemData.globEP = cellfun(@(E, P) E - P, problemData.globEP, globP, 'UniformOutput', false);

% Advection boundary edge integral (uu) without prescribed Dirichlet data for u (VI)
globPbdr = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~(problemData.g.markE0TbdrU | problemData.g.markE0TbdrRiemU | problemData.g.markE0TbdrCoupling), problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2});
globPRiem = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrRiemU, problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2});
problemData.globEP = cellfun(@(E, Pbdr, Priem) E - Pbdr - 0.5 * Priem, problemData.globEP, globPbdr, globPRiem, 'UniformOutput', false);

% Advection boundary edge integral (uu) with prescribed Dirichlet data for u (VI)
globJuu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
globJuuRiem = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiemU, @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = globJuu{1} + 0.5 * globJuuRiem{1};

% Advection boundary edge integral (uw) with prescribed Dirichlet data for u (VI)
problemData.globSbdr = assembleMatEdgeTetraPhiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.qOrd, problemData.hatQPerQuad);

% Advection boundary edge integral (gh) with prescribed Dirichlet data for h (VI)
globJh = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrH, @(x1,x2) hDCont(x1), problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
globJhRiem = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiemH, @(x1,x2) hDCont(x1), problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + problemData.gConst * (globJh{1} + 0.5 * globJhRiem{1});

% Diffusion element integral (IV)
problemData.globGR = assembleMatElemTetraDphiPhiFuncDisc(problemData.g, problemData.hatG, DDisc);

% Diffusion interior edge integral (V)
globR = assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);
problemData.globGR = cellfun(@(g,r) g - r, problemData.globGR, globR, 'UniformOutput', false);

% Diffusion boundary edge integral without prescribed Dirichlet data (V)
globRbdr = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrQ, problemData.hatRdiag, DDisc);
problemData.globGR = cellfun(@(g,r) g - r, problemData.globGR, globRbdr, 'UniformOutput', false);

% Diffusion boundary edge integral with prescribed Dirichlet data (V)
globJq = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrQ, qDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + globJq{1} + globJq{2};

%% Assembly of time-dependent vectors in flux and continuity equation.
% Boundary edge integrals with prescribed Dirichlet data for u (X, XII)
problemData.globJu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

%% Assembly of time-dependent matrices and vectors in continuity equation.
% Interior edge integrals (averaging uh/Hs) (XII)
problemData.tildeGlobP = problemData.fn_assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, problemData.g.markE0Tint, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);

% Boundary edge integrals with prescribed Dirichlet data for u or h and Riemann solver (XII)
problemData.tildeGlobPRiem = problemData.fn_assembleMatEdgeTetraPhiIntPhiIntFuncDisc1DIntNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, problemData.g.markE0TbdrRiemU | problemData.g.markE0TbdrRiemH, problemData.tildeHatPdiag);

% Boundary edge integrals with prescribed Dirichlet data and Riemann solver for u (XII)
problemData.globJuhRiem = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiemU & problemData.g.markE0TbdrRiemH, @(x1,x2) u1DCont(x1,x2) .* hDCont(x1), problemData.hSmoothV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Boundary edge integrals with prescribed Dirichlet data for w (XII)
problemData.globJw = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrW, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

%% Assembly of time-dependent matrices and vectors in free-surface equation.
% Element integrals (uh/Hs) in free surface equation (XIV)
problemData.barGlobG = problemData.fn_assembleMatElem1DDphiPhiFuncDiscHeight(barU1Disc, problemData.hSmoothQ0T1D, problemData.barHatG);

% Interior edge integrals (averaging uh/Hs) (XV)
problemData.barGlobP = problemData.fn_assembleMatV0T1DPhiPhiFuncDiscNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0Tint, problemData.barHatPdiag, problemData.barHatPoffdiag);

% Boundary edge integrals without prescribed Dirichlet data (XV)
markV0TbdrRiemUHUH = problemData.g.g1D.markV0TbdrRiemU | problemData.g.g1D.markV0TbdrRiemH | problemData.g.g1D.markV0TbdrRiemUH;
problemData.barGlobPbdr = problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0Tbdr & ~(problemData.g.g1D.markV0TbdrUH | markV0TbdrRiemUHUH), problemData.barHatPdiag);

% Boundary edge integrals with prescribed Dirichlet data and Riemann solver (XV)
problemData.barGlobPRiem = problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, markV0TbdrRiemUHUH, problemData.barHatPdiag);

% Boundary edge integrals with prescribed Dirichlet data for uh (XV)
problemData.barGlobJuh = zeros(problemData.g.g1D.numT * problemData.barN, 1);
problemData.barGlobJuhRiem = zeros(problemData.g.g1D.numT * problemData.barN, 1);
for n = 1 : 2
  markV0TbdrUH = problemData.g.g1D.markV0TbdrUH(:, n);
  markV0TbdrUHrep = logical(kron(markV0TbdrUH, true(problemData.barN, 1)));
  x1V0T = problemData.g.g1D.coordV0T(markV0TbdrUH, n, 1);
  problemData.barGlobJuh(markV0TbdrUHrep) = ( uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrUH, n) ) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n)';
  
  markV0TbdrRiemUH = problemData.g.g1D.markV0TbdrRiemUH(:, n);
  markV0TbdrRiemUHrep = logical(kron(markV0TbdrRiemUH, true(problemData.barN, 1)));
  x1V0T = problemData.g.g1D.coordV0T(markV0TbdrRiemUH, n, 1);
  problemData.barGlobJuhRiem(markV0TbdrRiemUHrep) = ( uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrRiemUH, n) ) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n)';
end % for n

%% Assembly of jump terms in Lax-Friedrichs Riemann-solver.
% 1D quadrature points
[Q,~] = quadRule1D(problemData.qOrd); 
numQuad1D = length(Q);

% Horizontal velocity from interior and neighboring element in quadrature points of edges
u1Q0E0Tint = { zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % cDisc{2} in quad points of edges
u1Q0E0TE0T = { zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % cDisc{2} of neighboring element in quad points of edges
u1Q0E0TbdrRiem = { zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % u1D in quad points of boundary edges with Riemann solver
for n = 1 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,n) * problemData.cDiscRK{nSubStep, 2}.', K * numQuad1D, 1);
  [Q1, Q2] = gammaMapTetra(n, Q);
  X1 = problemData.g.mapRef2Phy(1, Q1, Q2);
  X2 = problemData.g.mapRef2Phy(2, Q1, Q2);
  u1Q0E0TbdrRiem{n}(problemData.g.markE0TbdrRiem(:, n)) = u1DCont(X1(problemData.g.markE0TbdrRiem(:, n)), X2(problemData.g.markE0TbdrRiem(:, n)));
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,mapLocalEdgeTetra(n)) * problemData.cDiscRK{nSubStep, 2}.';
  u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', K * numQuad1D, 1);
end % for nn

% Jump in u in momentum equation (VI)
problemData.globKu = zeros(K * problemData.N, 1);
problemData.globKuRiem = zeros(K * problemData.N, 1);

% Jump in h in continuity equation (XII)
problemData.globKh = zeros(K * problemData.N, 1);
problemData.globKhRiem = zeros(K * problemData.N, 1);

% Jump in h in free-surface equation (XV)
problemData.barGlobKh = zeros(barK * problemData.barN, 1);
problemData.barGlobKhRiem = zeros(barK * problemData.barN, 1);

for n = 3 : 4
  nn1D = 5 - n; np1D = 5 - mapLocalEdgeTetra(n);
  
  [i, j] = find(problemData.g.g1D.markV0TbdrRiem(:, nn1D));
  hV0T1DbdrRiem = sparse(i, j, hDCont(problemData.g.g1D.coordV0T(problemData.g.g1D.markV0TbdrRiem(:, nn1D), nn1D, 1)), barK, 1);
  
  % Average and jump in height per edge
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,nn1D) + problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) + hV0T1DbdrRiem);
  hJmpE0T = problemData.g.g1D.markT2DT * ( ( problemData.hV0T1D(:,nn1D) - problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) - hV0T1DbdrRiem ) ./ problemData.hSmoothV0T1D(:,nn1D) );
  
  % Average and jump in horizontal velocity in quadrature points of edges
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n} + u1Q0E0TbdrRiem{n});
  u1JmpQ0E0T = u1Q0E0Tint{n} - u1Q0E0TE0T{n} - u1Q0E0TbdrRiem{n};
  
  % Eigenvalue of Jacobian of primitive numerical fluxes in quadrature points of edges
  lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  
  % Lax-Friedrichs Jump terms in quadrature points of edges
  hJmpLambdaQ0E0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D, 1));
  u1JmpLambdaQ0E0T = lambdaQ0E0T .* u1JmpQ0E0T;

  % Compute system vectors (VI, XII, XV)
  problemData.globKu = problemData.globKu + problemData.globS{n} * u1JmpLambdaQ0E0T;
  problemData.globKh = problemData.globKh + problemData.globS{n} * hJmpLambdaQ0E0T;
  problemData.barGlobKh = problemData.barGlobKh + problemData.barGlobS{n} * hJmpLambdaQ0E0T;
  
  problemData.globKuRiem = problemData.globKuRiem + problemData.globSuRiem{n} * u1JmpLambdaQ0E0T;
  problemData.globKhRiem = problemData.globKhRiem + problemData.globShRiem{n} * hJmpLambdaQ0E0T;
  problemData.barGlobKhRiem = problemData.barGlobKhRiem + problemData.barGlobSRiem{n} * hJmpLambdaQ0E0T;  
end % for n
end % function
