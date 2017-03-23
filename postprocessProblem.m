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
function problemData = postprocessProblem(problemData)
u1DCont = @(x1,x2) problemData.u1DCont(problemData.tEnd, x1, x2);
u2DCont = @(x1,x2) problemData.u2DCont(problemData.tEnd, x1, x2);
problemData.globJu = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrU, u1DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);
problemData.globJw = problemData.fn_assembleVecEdgeTrapPhiIntFuncContNu(problemData.g, problemData.g.markE0TbdrW, u2DCont, problemData.N, problemData.qOrd, problemData.basesOnQuad2D);

[Q,~] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);
heightV0T1D = problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,end), [4 3], 2) - problemData.g.coordV0T(problemData.g.g1D.idxT2D0T(:,1), [1 2], 2);

mapE0E = [2 1 4 3];
u1Q0E0Tint = cell(4,1); % cDisc{2} in quad points of edges
u1Q0E0TE0T = cell(4,1); % cDisc{2} of neighboring element in quad points of edges
for n = 1 : 4
  u1Q0E0Tint{n} = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,n) * problemData.cDiscRK{1, 2}.', problemData.g.numT * numQuad1D, 1);
  cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:,:,mapE0E(n)) * problemData.cDiscRK{1, 2}.';
  u1Q0E0TE0T{n} = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', problemData.g.numT * numQuad1D, 1);
end % for nn

problemData.globKh = zeros(problemData.g.numT * problemData.N, 1);
for n = 3 : 4
  hAvgE0T = 0.5 * problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,5-n) + problemData.g.g1D.markV0TV0T{5-n} * problemData.hV0T1D(:,5-mapE0E(n)) );
  hJmpE0T = problemData.g.g1D.markT2DT * ( problemData.hV0T1D(:,5-n) - problemData.g.g1D.markV0TV0T{5-n} * problemData.hV0T1D(:,5-mapE0E(n)) );
  u1AvgQ0E0T = 0.5 * (u1Q0E0Tint{n} + u1Q0E0TE0T{n});
  lambdaQ0E0T = 0.75 * abs(u1AvgQ0E0T) + 0.25 * sqrt( u1AvgQ0E0T .* u1AvgQ0E0T + 4 * problemData.gConst * kron(hAvgE0T, ones(numQuad1D,1)) );
  hJmpLambdaE0T = lambdaQ0E0T .* kron(hJmpE0T, ones(numQuad1D,1));
    
  problemData.globKh = problemData.globKh + problemData.globS{n} * hJmpLambdaE0T;
end % for n

problemData.tildeGlobP = assembleMatEdgeTrapPhiPhiFuncDisc1DNuHeight(problemData.g, problemData.g.g1D, problemData.cDiscRK{1, 1}, heightV0T1D, problemData.g.markE0Tint, problemData.tildeHatPdiag, problemData.tildeHatPoffdiag);

problemData.cDiscRK{1, 3} = reshape( (problemData.globHQup) \ ( problemData.globJu{1} + problemData.globJw{2} + problemData.globKh + ...
                                (problemData.globHQavg + problemData.tildeGlobP) * reshape(problemData.cDiscRK{1, 2}.', [], 1) ), ...
                             problemData.N, problemData.g.numT ).';

if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
  htEndCont = @(x1) problemData.hCont(problemData.tEnd, x1);
  u1tEndCont = @(x1,x2) problemData.u1Cont(problemData.tEnd, x1, x2);
  u2tEndCont = @(x1,x2) problemData.u2Cont(problemData.tEnd, x1, x2);

  problemData.error = [ computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{1, 1}, ...
                            htEndCont, problemData.qOrd + 1, problemData.basesOnQuad1D), ...
                        execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDiscRK{1, 2}, ...
                            u1tEndCont, problemData.qOrd + 1, problemData.basesOnQuad2D), ...
                        execin('darcyVert/computeL2ErrorTrap', problemData.g, problemData.cDiscRK{1, 3}, ...
                            u2tEndCont, problemData.qOrd + 1, problemData.basesOnQuad2D) ];

  fprintf('L2 errors of cDisc w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
end % if
if problemData.isVisGrid, execin('darcyVert/visualizeGridTrap', problemData.g); end
end % function

