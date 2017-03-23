% Preprocessing of the Runge-Kutta step.

%===============================================================================
%> @file advection/preprocessSubStep.m
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
DDisc = cellfun(@(c) problemData.fn_projectFuncCont2DataDiscTrap(problemData.g, @(x1,x2) c(t,x1,x2), problemData.N, problemData.qOrd, ...
        problemData.globM, problemData.basesOnQuad2D), problemData.DCont, 'UniformOutput', false);
problemData.globLu = reshape(problemData.fn_projectFuncCont2DataDiscTrap(problemData.g, @(x1,x2) problemData.fuCont(t,x1,x2), ...
                             problemData.N, problemData.qOrd, problemData.globM, problemData.basesOnQuad2D).', [], 1);
problemData.globLh = reshape(projectFuncCont2DataDisc1D(problemData.g.g1D, @(x1) problemData.fhCont(t,x1), ...
                             problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D).', [], 1);
%% Determine quadrature rule and mapping of local edge indices
[Q,~] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);
mapE0E = [2 1 4 3];
      
%% Create lookup tables for solution on quadrature points
              
u1Q0E0Tint = cell(4,1); % cDisc{2} in quad points of edges
u1Q0E0TE0T = cell(4,1); % cDisc{2} of neighboring element in quad points of edges
for n = 1 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,n) * problemData.cDiscRK{nSubStep, 2}.', problemData.g.numT * numQuad1D, 1);
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,mapE0E(n)) * problemData.cDiscRK{nSubStep, 2}.';
  u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', problemData.g.numT * numQuad1D, 1);
end % for nn

%% Compute depth averaged velocity
barU1Disc = { zeros(problemData.g.g1D.numT, problemData.barN), zeros(problemData.g.g1D.numT, problemData.barN) };
for s = 1 : 2
  for j = 1 : problemData.barN
    i = problemData.fn_mapTensorProductIndex(j, 1);
    barU1Disc{s}(:, j) = problemData.g.g1D.markT2DT.' * (problemData.cDiscRK{nSubStep, 2}(:,i) .* problemData.g.J0T{s}(:,2,2));
  end % for j
end % for s

%% Assembly of time-dependent global matrices.
% Advection and diffusion element integrals in momentum equation
problemData.globEP = problemData.fn_assembleMatElemTrapDphiPhiFuncDisc(problemData.g, problemData.hatG, problemData.cDiscRK{nSubStep, 2});%(nSubStep, 2:3));
% problemData.globEP = globE{1} + globE{2};
problemData.globGR = problemData.fn_assembleMatElemTrapDphiPhiFuncDisc(problemData.g, problemData.hatG, DDisc);

% Advection and diffusion interior edge integrals in momentum equation
globR = problemData.fn_assembleMatEdgeTrapPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, DDisc);
problemData.globGR = cellfun(@(g,r) g - r, problemData.globGR, globR, 'UniformOutput', false);
problemData.globEP = cellfun(@(e,p) e - p, problemData.globEP, problemData.fn_assembleMatEdgeTrapPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDiscRK{nSubStep, 2}), 'UniformOutput', false);
% globP = problemData.fn_assembleMatEdgeTrapPhiPhiFuncDiscNu(problemData.g, problemData.g.markE0Tint, problemData.hatRdiag, problemData.hatRoffdiag, problemData.cDiscRK(nSubStep, 2:3));
% problemData.globEP = problemData.globEP - globP{1} - globP{2};

% Jump term in Lax-Friedrichs solver
problemData.globKu = zeros(problemData.g.numT * problemData.N, 1);
problemData.globKh = zeros(problemData.g.numT * problemData.N, 1);
problemData.barGlobKh = zeros(problemData.g.g1D.numT * problemData.barN, 1);
for n = 3 : 4
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,5-n) + problemData.g.g1D.markV0TV0T{5-n} * problemData.hV0T1D(:,5-mapE0E(n)) );
  hJmpE0T = problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,5-n) - problemData.g.g1D.markV0TV0T{5-n} * problemData.hV0T1D(:,5-mapE0E(n)) );
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n});
  lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  hJmpLambdaE0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D,1));
    
  problemData.globKu = problemData.globKu + problemData.globS{n} * ( lambdaQ0E0T .* (u1Q0E0Tint{n} - u1Q0E0TE0T{n}) );
  problemData.globKh = problemData.globKh + problemData.globS{n} * hJmpLambdaE0T;
  problemData.barGlobKh = problemData.barGlobKh + (problemData.barGlobS{n} * hJmpLambdaE0T) ./ kron(problemData.heightV0T1D(:, 5-n), ones(problemData.barN, 1));
end % for n

% Interior edge flux in continuity and free surface equation
problemData.tildeGlobP = assembleMatEdgeTrapPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{nSubStep, 1}, problemData.heightV0T1D, problemData.g.markE0Tint, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);
problemData.barGlobP = assembleMatEdge1DPhiPhiFuncDiscNuHeight(problemData.g.g1D, barU1Disc, problemData.heightV0T1D, problemData.g.g1D.markV0Tint, problemData.barHatPdiag, problemData.barHatPoffdiag);

% Advection element integral in free surface equation
problemData.barGlobG = assembleMatElem1DDphiPhiFuncDiscHeight(barU1Disc, problemData.heightQ0T1D, problemData.barHatG);

%% Assembly of boundary contributions.
hDCont = @(x1) problemData.hDCont(t,x1);
u1DCont = @(x1,x2) problemData.u1DCont(t,x1,x2);
u2DCont = @(x1,x2) problemData.u2DCont(t,x1,x2);
qDCont = { @(x1,x2) problemData.q1DCont(t,x1,x2); @(x1,x2) problemData.q2DCont(t,x1,x2) };
uhDCont = @(x1) problemData.uhDCont(t, x1);

% Diffusion boundary terms in momentum equation
globRbdr = problemData.fn_assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrQ, problemData.hatRdiag, DDisc);
problemData.globGR = cellfun(@(g,r) g - r, problemData.globGR, globRbdr, 'UniformOutput', false);
globJtmp = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrQ, qDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = globJtmp{1} + globJtmp{2};

% Advection boundary terms in momentum equation
globJtmp = assembleVecEdgeTrapPhiIntFuncCont1DNu(problemData.g, problemData.g.markE0TbdrH, hDCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + problemData.gConst * globJtmp{1};
problemData.globEP = cellfun(@(e,p) e - p, problemData.globEP, problemData.fn_assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU, problemData.hatRdiag, problemData.cDiscRK{nSubStep, 2}), 'UniformOutput', false);
% globPbdr = problemData.fn_assembleMatEdgeTrapPhiIntPhiIntFuncDiscIntNu(problemData.g, problemData.g.markE0Tbdr & ~problemData.g.markE0TbdrU, problemData.hatRdiag, problemData.cDiscRK(nSubStep, 2:3));
% problemData.globEP = problemData.globEP - globPbdr{1} - globPbdr{2};

globJtmp = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, @(x1,x2) u1DCont(x1,x2).^2, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJ = problemData.globJ + globJtmp{1};
% globQPerQuad = assembleMatEdgeTrapPhiIntFuncDiscIntNuPerQuad(problemData.g, problemData.g.markE0TbdrU, problemData.cDiscRK{nSubStep, 3}, problemData.hatQPerQuad);
% for n = 1 : 4
%   [Q1, Q2] = execin('darcyVert/gammaMapTrap', n, Q);
%   u1DContQ0E0T = u1DCont(problemData.g.mapRef2Phy(1, Q1, Q2), problemData.g.mapRef2Phy(2, Q1, Q2));
%   problemData.globJ = problemData.globJ + globQPerQuad{n,2} * reshape(u1DContQ0E0T.', [], 1);
% end % for n
problemData.globSbdr = assembleMatEdgeTrapPhiIntPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.qOrd, problemData.hatQPerQuad);

% Boundary terms for horizontal velocity in continuity and flux equation
problemData.globJu = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Dirichlet boundary for vertical velocity in continuity equation
problemData.globJw = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrW, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

% Boundary terms in free surface equation
problemData.barGlobPbdr = assembleMatEdge1DPhiIntPhiIntFuncDiscIntNuHeight(problemData.g.g1D, barU1Disc, problemData.heightV0T1D, problemData.g.g1D.markV0Tbdr & ~problemData.g.g1D.markV0TbdrUH, problemData.barHatPdiag);
problemData.barGlobJuh = zeros(problemData.g.g1D.numT * problemData.barN, 1);
for n = 1 : 2
  markV0TbdrUH = problemData.g.g1D.markV0TbdrUH(:, n);
  markV0TbdrUHrep = logical(kron(markV0TbdrUH, true(problemData.barN, 1)));
  x1V0T = problemData.g.g1D.coordV0T(markV0TbdrUH, n, 1);
  problemData.barGlobJuh(markV0TbdrUHrep) = ( uhDCont(x1V0T) .* problemData.g.g1D.nuV0T(markV0TbdrUH, n) ) * problemData.basesOnQuad1D.phi0D{problemData.qOrd}(:, n)';
end % for n
end % function
