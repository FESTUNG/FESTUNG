% Performs all post-processing steps. Error evaluation for analytical problems.

%===============================================================================
%> @file sweVert/postprocessProblem.m
%>
%> @brief Performs all post-processing tasks. Error evaluation for analytical problems.
%===============================================================================
%>
%> @brief Performs all post-processing tasks. Error evaluation for analytical problems.
%>
%> This routine is called after the main loop.
%>
%> @param  problemData  A struct with problem parameters and solution
%>                      vectors. @f$[\text{struct}]@f$
%>
%> @retval problemData  A struct with all necessary parameters and definitions
%>                      for the problem description and precomputed fields.
%>                      @f$[\text{struct}]@f$
%>
%> This file is part of FESTUNG
%>
%> @copyright 2014-2017 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
function problemData = postprocessProblem(problemData)
%% Adapt mesh to final solution.
problemData = problemData.fn_adaptFreeSurface(problemData);

%% Compute vertical velocity.
u1DCont = @(x1,x2) problemData.u1DCont(problemData.tEnd, x1, x2);
u2DCont = @(x1,x2) problemData.u2DCont(problemData.tEnd, x1, x2);
problemData.globJu = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU & ~problemData.g.markE0TbdrCoupling, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJw = assembleVecEdgeTetraPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrW & ~problemData.g.markE0TbdrCoupling, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

[Q,~] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);
heightV0T1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), [4 3], 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), [1 2], 2);

u1Q0E0Tint = cell(4,1); % cDisc{2} in quad points of edges
u1Q0E0TE0T = cell(4,1); % cDisc{2} of neighboring element in quad points of edges
for n = 1 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,n) * problemData.cDiscRK{end, 2}.', problemData.g.numT * numQuad1D, 1);
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,mapLocalEdgeTetra(n)) * problemData.cDiscRK{end, 2}.';
  u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', problemData.g.numT * numQuad1D, 1);
end % for nn

problemData.globKh = zeros(problemData.g.numT * problemData.N, 1);
for n = 3 : 4
  nn1D = 5 - n; np1D = 5 - mapLocalEdgeTetra(n);
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,nn1D) + problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) );
  hJmpE0T = problemData.g.g1D.markT2DT * ( ( problemData.hV0T1D(:,nn1D) - problemData.g.g1D.markV0TV0T{nn1D} * problemData.hV0T1D(:,np1D) ) ./ problemData.hSmoothV0T1D(:,nn1D) );
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n});
  lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  hJmpLambdaE0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D,1));

  problemData.globKh = problemData.globKh + problemData.globS{n} * hJmpLambdaE0T;
end % for n

problemData.tildeGlobP = problemData.fn_assembleMatEdgeTetraPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{end, 1}, heightV0T1D, problemData.g.markE0Tint, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);

problemData.cDiscRK{end, 3} = reshape( (problemData.globHQup) \ ( problemData.globJu{1} + problemData.globJw{2} + problemData.globKh + ...
                                problemData.globJuCoupling{1} + problemData.globJwCoupling + ...
                                (problemData.globHQavg + problemData.tildeGlobP) * reshape(problemData.cDiscRK{end, 2}.', [], 1) ), ...
                             problemData.N, problemData.g.numT ).';

%% Evaluate errors (if analytical solution given).
if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
  htEndCont = @(x1) problemData.hCont(problemData.tEnd, x1);
  u1tEndCont = @(x1,x2) problemData.u1Cont(problemData.tEnd, x1, x2);
  u2tEndCont = @(x1,x2) problemData.u2Cont(problemData.tEnd, x1, x2);

  problemData.error = [ computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{end, 1}, ...
                            htEndCont, problemData.qOrd + 1, problemData.basesOnQuad1D), ...
                        computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 2}, ...
                            u1tEndCont, problemData.qOrd + 1, problemData.basesOnQuad2D), ...
                        computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 3}, ...
                            u2tEndCont, problemData.qOrd + 1, problemData.basesOnQuad2D) ];

  fprintf('L2 errors of cDisc w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
end % if

%% Visualize final state.
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

if problemData.isVisSol
  cLagr = { projectDataDisc2DataLagrTensorProduct(problemData.cDiscRK{end, 2}), ...
            projectDataDisc2DataLagrTensorProduct(problemData.cDiscRK{end, 3}) };
  visualizeDataLagrTetra(problemData.g, cLagr, {'u1', 'u2'}, problemData.outputBasename, problemData.numSteps, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
end % function

