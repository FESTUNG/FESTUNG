% Preprocessing of the Runge-Kutta step.

%===============================================================================
%> @file
%>
%> @brief Preprocessing of the Runge-Kutta step.
%===============================================================================
%>
%> @brief Preprocessing of the Runge-Kutta step.
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
%> This routine is executed first in each loop iteration.
%> It takes care of the assembly of time-dependent matrices and right hand
%> side vectors. Furthermore, it evaluates the boundary conditions for the 
%> current Runge-Kutta step.
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
%> @retval problemData  The input struct enriched with preprocessed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
function problemData = preprocessSubStep(problemData, nStep, nSubStep) %#ok<INUSL>
t = problemData.t(nSubStep);
K = problemData.g.numT;
barK = problemData.g.g1D.numT;

[Q,~] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);

%% Algebraic right hand sides and boundary conditions.
fuCont = problemData.fuCont;
fhCont = problemData.fhCont;

hDCont = problemData.hDCont;
u1DCont = problemData.u1DCont;
u2DCont = problemData.u2DCont;
qDCont = { @(x1,x2) problemData.q1DCont(t,x1,x2); @(x1,x2) problemData.q2DCont(t,x1,x2) };

%% L2-projections of algebraic coefficients and right hand side.
if iscell(problemData.DCont)
  DDisc = cellfun(@(c) projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) c(t,x1,x2), problemData.qOrd, ...
                         problemData.globM, problemData.basesOnQuad2D), problemData.DCont, 'UniformOutput', false);
else
  DDisc = projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) problemData.DCont(t,x1,x2), problemData.qOrd, ...
                         problemData.globM, problemData.basesOnQuad2D);
end % if
problemData.globLu = reshape(projectFuncCont2DataDiscTetra(problemData.g, @(x1,x2) fuCont(t,x1,x2), problemData.qOrd, problemData.globM, problemData.basesOnQuad2D).', [], 1);
problemData.globLh = reshape(projectFuncCont2DataDisc1D(problemData.g.g1D, @(x1) fhCont(t,x1), problemData.qOrd, problemData.hatBarM, problemData.basesOnQuad1D).', [], 1);

%% Compute depth integrated velocity and water height in surface nodes.
barU1Disc = { zeros(problemData.g.g1D.numT, problemData.barN), zeros(problemData.g.g1D.numT, problemData.barN) };
for s = 1 : 2
  for j = 1 : problemData.barN
    i = mapTensorProductIndex(j, 1);
    barU1Disc{s}(:, j) = problemData.g.g1D.markT2DT.' * (problemData.cDiscRK{nSubStep, 2}(:,i) .* problemData.g.J0T{s}(:,2,2));
  end % for j
end % for s
hSmooth_hV0T1D = problemData.hSmoothV0T1D ./ (problemData.cDiscRK{nSubStep, 1} * problemData.basesOnQuad1D.phi0D{problemData.qOrd});
hSmooth_hDV0T1D = problemData.hSmoothV0T1D ./ hDCont(t, problemData.g.g1D.coordV0T(:,:,1));

%% Assembly of time-dependent matrices and vectors in momentum equation.
% Advection element integral (II)
globE = assembleMatElemDphiPhiFuncDisc(problemData.g, problemData.hatG, problemData.cDiscRK{nSubStep, 2});

% Advection interior edge integrals (averaging uu, uw) (VI)
markE0T = problemData.g.markE0Tint | (problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU);
globP = assembleMatEdgePhiPhiFuncDiscNu(problemData.g, markE0T, problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDiscRK{nSubStep, 2});

% Advection boundary edge integral (uu) without prescribed Dirichlet data for u (VI)
% TEST bottom friction: markE0T = problemData.g.markE0TbdrTop | (problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU) | (~problemData.isCoupling & problemData.isBottomFriction & problemData.g.markE0TbdrBot);
markE0T = problemData.g.markE0TbdrTop | (problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU);
globPbdr = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(problemData.g, markE0T, problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2});

% Advection boundary edge integral (uu) with prescribed Dirichlet data for u (VI)
% TEST bottom friction: globJbot = assembleVecEdgePhiIntFuncContNu(problemData.g, ~problemData.isCoupling & ~problemData.isBottomFriction & problemData.g.markE0TbdrBot, { @(x1,x2) u1DCont(t,x1,x2).^2, @(x1,x2) u1DCont(t,x1,x2) .* u2DCont(t,x1,x2) }, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
globJbot = assembleVecEdgePhiIntFuncContNu(problemData.g, ~problemData.isCoupling & problemData.g.markE0TbdrBot, { @(x1,x2) u1DCont(t,x1,x2).^2, @(x1,x2) u1DCont(t,x1,x2) .* u2DCont(t,x1,x2) }, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
globJuu = assembleVecEdgePhiIntFuncContNu(problemData.g, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU, @(x1,x2) u1DCont(t,x1,x2).^2, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
globJuuRiem = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU, @(x1,x2) u1DCont(t,x1,x2).^2, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);

% Advection boundary edge integral (gh) with prescribed Dirichlet data for h (VI)
globJh = assembleVecEdgePhiIntFuncContNu(problemData.g, ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrH, @(x1,x2) hDCont(t, x1), problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
globJhRiem = assembleVecEdgePhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrH, @(x1,x2) hDCont(t, x1), problemData.N, problemData.basesOnQuad2D, problemData.qOrd);

% Diffusion element integral (IV)
globG = assembleMatElemDphiPhiFuncDisc(problemData.g, problemData.hatG, DDisc);

% Diffusion interior edge integral (V)
globR = assembleMatEdgePhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);

% Diffusion boundary edge integral without prescribed Dirichlet data (V)
markE0T = (~problemData.isBottomFriction & problemData.g.markE0TbdrBot) | (problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrQ);
globRbdr = assembleMatEdgePhiIntPhiIntFuncDiscIntNu(problemData.g, markE0T, problemData.hatRdiag, DDisc);

% Diffusion boundary edge integral with prescribed Dirichlet data (V)
markE0T = problemData.g.markE0TbdrQ | problemData.g.markE0TbdrTop;
globJq = assembleVecEdgePhiIntFuncContNu(problemData.g, markE0T, qDCont, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);

%% Assembly of bottom friction term
globJfric = zeros(K * problemData.N, 1);
if problemData.isBottomFriction
  u1Q0E0T = zeros(K, 4, numQuad1D);
  u1Q0E0T(problemData.g.markE0TbdrBot, 1, :) = problemData.cDiscRK{nSubStep, 2}(problemData.g.markE0TbdrBot, :) * problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 1).';
   
  switch problemData.bottomFriction
    case 'linear'
      valOnQuad = -problemData.CfConst * u1Q0E0T;
      
    case 'quadratic'
      valOnQuad = -problemData.CfConst * abs(u1Q0E0T) .* u1Q0E0T;
      
    otherwise
      error('Invalid bottom friction type')
  end % switch
  
  globJfric = assembleVecEdgePhiIntVal(problemData.g, problemData.g.markE0TbdrBot, valOnQuad, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
end % if

%% Put together matrices and vectors
problemData.globEP = cellfun(@(E, P, Pbdr) E - P - Pbdr, globE, globP, globPbdr, 'UniformOutput', false);
problemData.globGR = cellfun(@(G, R, Rbdr) G - R - Rbdr, globG, globR, globRbdr, 'UniformOutput', false);
problemData.globJ = globJfric + globJbot{1} + globJbot{2} + globJuu{1} + 0.5 * globJuuRiem{1} + problemData.gConst * (globJh{1} + 0.5 * globJhRiem{1}) + globJq{1} + globJq{2};

%% Assembly of time-dependent vectors in flux and continuity equation.
% Boundary edge integrals with prescribed Dirichlet data for u (X, XII)
% TEST bottom friction: markE0T = (~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU) | (~problemData.isCoupling & ~problemData.isBottomFriction & problemData.g.markE0TbdrBot);
markE0T = (~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU) | (~problemData.isCoupling & problemData.g.markE0TbdrBot);
problemData.globJu = assembleVecEdgePhiIntFuncContNu(problemData.g, markE0T, { @(x1,x2) u1DCont(t,x1,x2), @(x1,x2) u2DCont(t,x1,x2) }, problemData.N, problemData.basesOnQuad2D, problemData.qOrd);
% TEST bottom friction: markE0T = problemData.g.markE0TbdrU | (~problemData.isCoupling & ~problemData.isBottomFriction & problemData.g.markE0TbdrBot);
markE0T = problemData.g.markE0TbdrU | (~problemData.isCoupling & problemData.g.markE0TbdrBot);
problemData.globJuFlux = assembleVecEdgePhiIntFuncContNu(problemData.g, markE0T, @(x1,x2) u1DCont(t,x1,x2), problemData.N, problemData.basesOnQuad2D, problemData.qOrd);

%% Assembly of time-dependent matrices and vectors in continuity equation.
% Interior edge integrals (averaging uh/Hs) (XII)
markE0T = (problemData.g.markE0Tint & problemData.g.markE0Tv) | (problemData.g.markE0TbdrRiem & (problemData.g.markE0TbdrU | problemData.g.markE0TbdrH));
problemData.globVeeP = problemData.fn_assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, markE0T, problemData.hatVeePdiag, problemData.hatVeePoffdiag);

% Boundary edge integrals with no prescribed Dirichlet and Riemann solver for u and h
markE0T = problemData.g.markE0Tbdr & ~(problemData.g.markE0TbdrU | problemData.g.markE0TbdrH);
problemData.globVeePbdr = problemData.fn_assembleMatEdgeTetraPhiIntPhiIntFuncDisc1DIntNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.hSmoothV0T1D, markE0T, problemData.hatVeePdiag);

% Boundary edge integrals with prescribed Dirichlet data and Riemann solver for u and h (XII)
problemData.globJuhRiem = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & problemData.g.markE0TbdrH, @(x1,x2) u1DCont(t,x1,x2) .* hDCont(t,x1), problemData.hSmoothV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Boundary edge integrals with prescribed Dirichlet Data and no Riemann solver for u or h but not both (XII)
markE0T = problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH;
problemData.globVeeJu = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, markE0T, @(x1,x2) u1DCont(t,x1,x2), hSmooth_hV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
markE0T = problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrRiem & ~problemData.g.markE0TbdrU & problemData.g.markE0TbdrH;
problemData.globVeeJh = problemData.fn_assembleVecEdgeTetraPhiIntFuncDiscIntHeightNu(problemData.g, problemData.g.g1D, markE0T, problemData.cDiscRK{nSubStep, 2}, hSmooth_hDV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Boundary edge integrals with prescribed Dirichlet Data and Riemann solver for u or h but not both (XII)
markE0T = problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH;
problemData.globVeeJuRiem = problemData.fn_assembleVecEdgeTetraPhiIntFuncContHeightNu(problemData.g, problemData.g.g1D, markE0T, @(x1,x2) u1DCont(t,x1,x2), hSmooth_hV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
markE0T = problemData.g.markE0TbdrRiem & ~problemData.g.markE0TbdrU & problemData.g.markE0TbdrH;
problemData.globVeeJhRiem = problemData.fn_assembleVecEdgeTetraPhiIntFuncDiscIntHeightNu(problemData.g, problemData.g.g1D, markE0T, problemData.cDiscRK{nSubStep, 2}, hSmooth_hDV0T1D, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

%% Assembly of time-dependent matrices and vectors in free-surface equation.
% Element integrals (uh/Hs) in free surface equation (XIV)
problemData.globBarG = problemData.fn_assembleMatElem1DDphiPhiFuncDiscHeight(barU1Disc, problemData.hSmoothQ0T1D, problemData.hatBarG);

% Interior edge integrals (averaging uh/Hs) (XV)
markV0TbdrUoH = problemData.g.g1D.markV0TbdrU | problemData.g.g1D.markV0TbdrH;
markV0T = problemData.g.g1D.markV0Tint | (problemData.g.g1D.markV0TbdrRiem & markV0TbdrUoH);
problemData.globBarP = problemData.fn_assembleMatV0T1DPhiPhiFuncDiscNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, markV0T, problemData.hatBarPdiag, problemData.hatBarPoffdiag);

% Boundary edge integrals without prescribed Dirichlet data (XV)
markV0T = problemData.g.g1D.markV0Tbdr & ~markV0TbdrUoH;
problemData.globBarPbdr = problemData.fn_assembleMatV0T1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.hSmoothV0T1D, markV0T, problemData.hatBarPdiag);

% Boundary edge integrals with prescribed Dirichlet data for uh (XV)
markV0T = ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & problemData.g.markE0TbdrH;
problemData.globBarJuh = problemData.fn_assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, markV0T, @(x1,x2) u1DCont(t,x1,x2) .* hDCont(t,x1), problemData.hSmoothV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
markV0T = problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & problemData.g.markE0TbdrH;
problemData.globBarJuhRiem = problemData.fn_assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, markV0T, @(x1,x2) u1DCont(t,x1,x2) .* hDCont(t,x1), problemData.hSmoothV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

markV0T = ~problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH;
problemData.globBarJu = problemData.fn_assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, markV0T, @(x1,x2) u1DCont(t,x1,x2), hSmooth_hV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
markV0T = problemData.g.markE0TbdrRiem & problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrH;
problemData.globBarJuRiem = problemData.fn_assembleVecV0T1DPhiIntFuncContNuHeight(problemData.g, problemData.g.g1D, markV0T, @(x1,x2) u1DCont(t,x1,x2), hSmooth_hV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

markV0T = ~problemData.g.g1D.markV0TbdrRiem & ~problemData.g.g1D.markV0TbdrU & problemData.g.g1D.markV0TbdrH;
problemData.globBarJh = problemData.fn_assembleVecV0T1DPhiIntFuncDiscIntNuHeight(problemData.g.g1D, markV0T, barU1Disc, hSmooth_hDV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);
markV0T = problemData.g.g1D.markV0TbdrRiem & ~problemData.g.g1D.markV0TbdrU & problemData.g.g1D.markV0TbdrH;
problemData.globBarJhRiem = problemData.fn_assembleVecV0T1DPhiIntFuncDiscIntNuHeight(problemData.g.g1D, markV0T, barU1Disc, hSmooth_hDV0T1D, problemData.barN, problemData.qOrd, problemData.basesOnQuad1D);

%% Assembly of jump terms in Lax-Friedrichs Riemann-solver.

% Horizontal velocity from interior and neighboring element in quadrature points of edges
u1Q0E0Tint = zeros(K * numQuad1D, 2); % cDisc{2} in quad points of edges
u1Q0E0TE0T = zeros(K * numQuad1D, 2); % cDisc{2} of neighboring element in quad points of edges
u1Q0E0TbdrRiem = zeros(K * numQuad1D, 2); % u1D in quad points of boundary edges with Riemann solver
for n = 3 : 4
  u1Q0E0Tint(:,n-2) = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, n) * problemData.cDiscRK{nSubStep, 2}.', K * numQuad1D, 1);
  
  markQ0E0TbdrRiem = logical(kron(problemData.g.markE0TbdrRiem(:, n) & problemData.g.markE0TbdrU(:, n), true(numQuad1D, 1)));
  [Q1, Q2] = gammaMapTetra(n, Q);
  X1 = reshape(problemData.g.mapRef2Phy(1, Q1, Q2).', K * numQuad1D, 1);
  X2 = reshape(problemData.g.mapRef2Phy(2, Q1, Q2).', K * numQuad1D, 1);
  u1Q0E0TbdrRiem(markQ0E0TbdrRiem, n-2) = u1DCont(t, X1(markQ0E0TbdrRiem), X2(markQ0E0TbdrRiem));
  
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, mapLocalEdgeTetra(n)) * problemData.cDiscRK{nSubStep, 2}.';
  u1Q0E0TE0T(:, n-2) = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', K * numQuad1D, 1);
end % for nn

% Jump in u in momentum equation (VI)
problemData.globKu = zeros(K * problemData.N, 1);

% Jump in h in continuity equation (XII)
problemData.globKh = zeros(K * problemData.N, 1);

% Jump in h in free-surface equation (XV)
problemData.globBarKh = zeros(barK * problemData.barN, 1);

for n = 3 : 4
  nn1D = 5 - n; np1D = 5 - mapLocalEdgeTetra(n);
  markV0TbdrRiem = problemData.g.g1D.markV0TbdrRiem(:, nn1D) & problemData.g.g1D.markV0TbdrH(:, nn1D);
  
  [i, j] = find(markV0TbdrRiem);
  hV0T1DbdrRiem = sparse(i, j, hDCont(t, problemData.g.g1D.coordV0T(markV0TbdrRiem, nn1D, 1)), barK, 1);
  
  % Average and jump in height per edge
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,nn1D) + problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) + hV0T1DbdrRiem);
  hJmpE0T = problemData.g.g1D.markT2DT * ( ( problemData.hV0T1D(:,nn1D) - problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) - hV0T1DbdrRiem ) ./ problemData.hSmoothV0T1D(:,nn1D) );
  
  % Average and jump in horizontal velocity in quadrature points of edges
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint(:, n-2) + u1Q0E0TE0T(:, n-2) + u1Q0E0TbdrRiem(:, n-2));
  u1JmpQ0E0T = u1Q0E0Tint(:, n-2) - u1Q0E0TE0T(:, n-2) - u1Q0E0TbdrRiem(:, n-2);
  
  % Eigenvalue of Jacobian of primitive numerical fluxes in quadrature points of edges
  lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  
  % Lax-Friedrichs Jump terms in quadrature points of edges
  hJmpLambdaQ0E0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D, 1));
  u1JmpLambdaQ0E0T = lambdaQ0E0T .* u1JmpQ0E0T;

  % Compute system vectors (VI, XII, XV)
  problemData.globKu = problemData.globKu + problemData.globSu{n} * u1JmpLambdaQ0E0T;
  problemData.globKh = problemData.globKh + problemData.globSh{n} * hJmpLambdaQ0E0T;
  problemData.globBarKh = problemData.globBarKh + problemData.globBarS{n} * hJmpLambdaQ0E0T; 
end % for n
end % function
