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

%% L2-projections of algebraic coefficients and right hand side.
DDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
                       problemData.globM, problemData.basesOnQuad2D), problemData.DCont, 'UniformOutput', false);
problemData.globLu = reshape(projectFuncCont2DataDiscTetra(problemData.g, fuCont, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D).', [], 1);
problemData.globLh = reshape(projectFuncCont2DataDisc1D(problemData.g.g1D, fhCont, problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D).', [], 1);

%% Compute depth integrated velocity and water height in surface nodes.
barU1Disc = { zeros(problemData.g.g1D.numT, problemData.barN), zeros(problemData.g.g1D.numT, problemData.barN) };
for s = 1 : 2
  for j = 1 : problemData.barN
    i = mapTensorProductIndex(j, 1);
    barU1Disc{s}(:, j) = problemData.g.g1D.markT2DT.' * (problemData.cDiscRK{nSubStep, 2}(:,i) .* problemData.g.J0T{s}(:,2,2));
  end % for j
end % for s
hSmooth_hV0T1D = problemData.hSmoothV0T1D ./ (problemData.cDiscRK{nSubStep, 1} * problemData.basesOnQuad1D.phi0D{problemData.qOrd});
hSmooth_hDV0T1D = problemData.hSmoothV0T1D ./ hDCont(problemData.g.g1D.coordV0T(:,:,1));

%% Assembly of time-dependent matrices and vectors in momentum equation.
% Advection element integral (II)
globE = assembleMatElemTetraDphiPhiFuncDisc(problemData.g, problemData.hatG, problemData.cDiscRK{nSubStep, 2});

% Advection interior edge integrals (averaging uu, uw) (VI)
globP = assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDiscRK{nSubStep, 2});

% Advection boundary edge integral (uu) without prescribed Dirichlet data for u (VI)
globPtop = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrTop, problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2});
globPbdr = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU, problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2});
globPRiem = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU, problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2});

% Advection boundary edge integral (uu) with prescribed Dirichlet data for u (VI)
globJbot = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrBot, { @(x1,x2) u1DCont(x1,x2).^2, @(x1,x2) u1DCont(x1,x2) .* u2DCont(x1,x2) }, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
globJuu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU, @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
globJuuRiem = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU, @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Advection boundary edge integral (uw) with prescribed Dirichlet data for u (VI)
% problemData.globSbdr = assembleMatEdgeTetraPhiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.qOrd, problemData.hatQPerQuad);
% problemData.globSriem = assembleMatEdgeTetraPhiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiemU, u1DCont, problemData.qOrd, problemData.hatQPerQuad);

% Advection boundary edge integral (gh) with prescribed Dirichlet data for h (VI)
globJh = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrH, @(x1,x2) hDCont(x1), problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
globJhRiem = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrH, @(x1,x2) hDCont(x1), problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Diffusion element integral (IV)
globG = assembleMatElemTetraDphiPhiFuncDisc(problemData.g, problemData.hatG, DDisc);

% Diffusion interior edge integral (V)
globR = assembleMatEdgeTetraPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);

% Diffusion boundary edge integral without prescribed Dirichlet data (V)
globRbot = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0TbdrBot, problemData.hatRdiag, DDisc);
globRbdr = assembleMatEdgeTetraPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrQ, problemData.hatRdiag, DDisc);

% Diffusion boundary edge integral with prescribed Dirichlet data (V)
globJq = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrQ, qDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Put together matrices
problemData.globEP = cellfun(@(E, P, Ptop, Pbdr, PRiem) E - P - Ptop - Pbdr - 0.5 * PRiem, globE, globP, globPtop, globPbdr, globPRiem, 'UniformOutput', false);
problemData.globGR = cellfun(@(G, R, Rbot, Rbdr) G - R - Rbot - Rbdr, globG, globR, globRbot, globRbdr, 'UniformOutput', false);
problemData.globJ = globJbot{1} + globJbot{2} + globJuu{1} + 0.5 * globJuuRiem{1} + problemData.gConst * (globJh{1} + 0.5 * globJhRiem{1}) + globJq{1} + globJq{2};

%% Assembly of time-dependent vectors in flux and continuity equation.
% Boundary edge integrals with prescribed Dirichlet data for u (X, XII)
problemData.globJu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJuFlux = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJuBot = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrBot, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

%% Assembly of time-dependent matrices and vectors in continuity equation.
% Interior edge integrals (averaging uh/Hs) (XII)
problemData.tildeGlobP = problemData.fn_assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, problemData.g.markE0Tint & problemData.g.markE0Tv, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);

% Boundary edge integrals with prescribed Dirichlet data for u or h and Riemann solver (XII)
problemData.tildeGlobPRiem = problemData.fn_assembleMatEdgeTetraPhiIntPhiIntFuncDisc1DIntNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, problemData.g.markE0TbdrRiem & (problemData.g.markE0TbdrU | problemData.g.markE0TbdrH), problemData.tildeHatPdiag);

% Boundary edge integrals with prescribed Dirichlet data and Riemann solver for u and h (XII)
problemData.globJuhRiem = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & problemData.g.markE0TbdrH, @(x1,x2) u1DCont(x1,x2) .* hDCont(x1), problemData.hSmoothV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Boundary edge integrals with prescribed Dirichlet data for w (XII)
problemData.globJw = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrBot, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Boundary edge integrals with prescribed Dirichlet Data and Riemann solver for u or h but not both (XII)
problemData.tildeGlobJuRiem = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH, u1DCont, hSmooth_hV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.tildeGlobJhRiem = assembleVecEdgeTetraPhiIntFuncDiscIntHeightNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiem & ~problemData.g.markE0TbdrU & problemData.g.markE0TbdrH, problemData.cDiscRK{nSubStep, 2}, hSmooth_hDV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

%% Assembly of time-dependent matrices and vectors in free-surface equation.
% Element integrals (uh/Hs) in free surface equation (XIV)
problemData.barGlobG = problemData.fn_assembleMatElem1DDphiPhiFuncDiscHeight(barU1Disc, problemData.hSmoothQ0T1D, problemData.barHatG);

% Interior edge integrals (averaging uh/Hs) (XV)
problemData.barGlobP = problemData.fn_assembleMatV0T1DPhiPhiFuncDiscNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0Tint, problemData.barHatPdiag, problemData.barHatPoffdiag);

% Boundary edge integrals without prescribed Dirichlet data (XV)
markV0TbdrUoH = problemData.g.g1D.markV0TbdrU | problemData.g.g1D.markV0TbdrH;
problemData.barGlobPbdr = problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0Tbdr & ~markV0TbdrUoH, problemData.barHatPdiag);

% Boundary edge integrals with prescribed Dirichlet data and Riemann solver (XV)
problemData.barGlobPRiem = problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, problemData.g.g1D.markV0TbdrRiem & markV0TbdrUoH, problemData.barHatPdiag);

% Boundary edge integrals with prescribed Dirichlet data for uh (XV)
problemData.barGlobJuh = assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & problemData.g.markE0TbdrH, @(x1,x2) u1DCont(x1,x2) .* hDCont(x1), problemData.hSmoothV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barGlobJuhRiem = assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & problemData.g.markE0TbdrH, @(x1,x2) u1DCont(x1,x2) .* hDCont(x1), problemData.hSmoothV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

problemData.barGlobJu = assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH, u1DCont, hSmooth_hV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barGlobJuRiem = assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH, u1DCont, hSmooth_hV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

% barGlobJh = zeros(problemData.g.g1D.numT, problemData.barN);
% [~, W] = quadRule1D(problemData.qOrd);
% markE0T = ~problemData.g.markE0TbdrRiem & ~problemData.g.markE0TbdrU & problemData.g.markE0TbdrH;
% for n = 3 : 4
%   areaNuE0T = markE0T(:, n) .* problemData.g.areaE0T(:, n) .* problemData.g.nuE0T(:, n, 1);
%   funcIntegrated = (problemData.cDiscRK{nSubStep, 2} * problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, n).') * W(:);
%   barGlobJh = barGlobJh + ((problemData.g.g1D.markT2DT.' * (areaNuE0T .* funcIntegrated)) ./ hSmooth_hDV0T1D(:, 5-n)) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, 5-n).';
% end % for
% barGlobJh = reshape(barGlobJh.', [], 1);

problemData.barGlobJh = assembleVecV0T1DPhiIntFuncDiscIntNuHeight(problemData.g.g1D, ~problemData.g.g1D.markV0TbdrRiem & ~problemData.g.g1D.markV0TbdrU & problemData.g.g1D.markV0TbdrH, barU1Disc, hSmooth_hDV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
problemData.barGlobJhRiem = assembleVecV0T1DPhiIntFuncDiscIntNuHeight(problemData.g.g1D, problemData.g.g1D.markV0TbdrRiem & ~problemData.g.g1D.markV0TbdrU & problemData.g.g1D.markV0TbdrH, barU1Disc, hSmooth_hDV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
% barGlobJuh = zeros(problemData.g.g1D.numT * problemData.barN, 1);
% barGlobJuhRiem = zeros(problemData.g.g1D.numT * problemData.barN, 1);
% for n = 1 : 2
%   markV0TbdrUH = problemData.g.g1D.markV0TbdrUH(:, n);
%   markV0TbdrUHrep = logical(kron(markV0TbdrUH, true(problemData.barN, 1)));
%   x1V0T = problemData.g.g1D.coordV0T(markV0TbdrUH, n, 1);
%   barGlobJuh(markV0TbdrUHrep) = barGlobJuh(markV0TbdrUHrep) + ( uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrUH, n) ) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n)';
%   
%   markV0TbdrRiemUH = problemData.g.g1D.markV0TbdrRiemUH(:, n);
%   markV0TbdrRiemUHrep = logical(kron(markV0TbdrRiemUH, true(problemData.barN, 1)));
%   x1V0T = problemData.g.g1D.coordV0T(markV0TbdrRiemUH, n, 1);
%   barGlobJuhRiem(markV0TbdrRiemUHrep) = barGlobJuhRiem(markV0TbdrRiemUHrep) + ( (uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrRiemUH, n)) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n).' ).';
% end % for n

%% Assembly of jump terms in Lax-Friedrichs Riemann-solver.
% 1D quadrature points
[Q,~] = quadRule1D(problemData.qOrd); 
numQuad1D = length(Q);

% Horizontal velocity from interior and neighboring element in quadrature points of edges
u1Q0E0Tint = { []; []; zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % cDisc{2} in quad points of edges
u1Q0E0TE0T = { []; []; zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % cDisc{2} of neighboring element in quad points of edges
u1Q0E0TbdrRiem = { []; []; zeros(K * numQuad1D, 1); zeros(K * numQuad1D, 1) }; % u1D in quad points of boundary edges with Riemann solver
for n = 3 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, n) * problemData.cDiscRK{nSubStep, 2}.', K * numQuad1D, 1);
  
  markQ0E0TbdrRiem = logical(kron(problemData.g.markE0TbdrRiem(:, n) & problemData.g.markE0TbdrU(:, n), true(numQuad1D, 1)));
  [Q1, Q2] = gammaMapTetra(n, Q);
  X1 = reshape(problemData.g.mapRef2Phy(1, Q1, Q2).', K * numQuad1D, 1);
  X2 = reshape(problemData.g.mapRef2Phy(2, Q1, Q2).', K * numQuad1D, 1);
  u1Q0E0TbdrRiem{n}(markQ0E0TbdrRiem) = u1DCont(X1(markQ0E0TbdrRiem), X2(markQ0E0TbdrRiem));
  
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, mapLocalEdgeTetra(n)) * problemData.cDiscRK{nSubStep, 2}.';
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
  markV0TbdrRiem = problemData.g.g1D.markV0TbdrRiem(:, nn1D) & problemData.g.g1D.markV0TbdrH(:, nn1D);
  
  [i, j] = find(markV0TbdrRiem);
  hV0T1DbdrRiem = sparse(i, j, hDCont(problemData.g.g1D.coordV0T(markV0TbdrRiem, nn1D, 1)), barK, 1);
  
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
