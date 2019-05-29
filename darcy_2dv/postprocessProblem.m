% Performs all post-processing steps. Error evaluation for analytical problems.

%===============================================================================
%> @file
%>
%> @brief Performs all post-processing tasks. Error evaluation for analytical problems.
%===============================================================================
%>
%> @brief Performs all post-processing tasks. Error evaluation for analytical problems.
%>
%> This routine is called after the main loop.
%>
%> If analytical solutions are availabe, it evaluates the L2-errors.
%> Furthermore, the final state is written to file, allowing to restart the 
%> computation from the final state.
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
function problemData = postprocessProblem(problemData)
K = problemData.g.numT;
N = problemData.N;
  
hDisc = reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)';
q1Disc = reshape(problemData.sysY(1 : K*N), N, K)';
q2Disc = reshape(problemData.sysY(K*N+1 : 2*K*N), N, K)';

if all(isfield(problemData, { 'hCont', 'q1Cont', 'q2Cont' }))
  hCont = @(x1,x2) problemData.hCont(problemData.tEnd, x1, x2);
  q1Cont = @(x1,x2) problemData.q1Cont(problemData.tEnd, x1, x2);
  q2Cont = @(x1,x2) problemData.q2Cont(problemData.tEnd, x1, x2);

  problemData.error = [ computeL2ErrorQuadri(problemData.g, hDisc, hCont, ...
                          problemData.qOrd+1, problemData.basesOnQuad), ...
                        computeL2ErrorQuadri(problemData.g, q1Disc, q1Cont, ...
                          problemData.qOrd+1, problemData.basesOnQuad), ...
                        computeL2ErrorQuadri(problemData.g, q2Disc, q2Cont, ...
                          problemData.qOrd+1, problemData.basesOnQuad) ];

  fprintf('L2 errors of h, q1, q2 w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
end % if

%% Norms for stability estimate
% [Q, W] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);

% hDisc = problemData.sysY(2*K*N+1 : 3*K*N);
% normH = hDisc.' * (problemData.globM * hDisc);

% q1Disc = problemData.sysY(1 : K*N);
% normQ1 = q1Disc.' * (problemData.globM * q1Disc);

% q2Disc = problemData.sysY(K*N+1 : 2*K*N);
% normQ2 = q2Disc.' * (problemData.globM * q2Disc);

% normJumpH = 0;
% hDisc = reshape(hDisc, N, K)';
% for n = 2 : 3
%   hQ0E0Tint = reshape(problemData.basesOnQuad.phi1D{problemData.qOrd}(:, :, n) * bsxfun( ...
%                       @times, problemData.g.markE0Tint(:, n), hDisc).', ...
%                   K * numQuad1D, 1);
%   cDiscThetaPhi = problemData.basesOnQuad.phi1D{problemData.qOrd}(:, :, mapLocalEdgeIndexQuadri(3)) * hDisc.';
%   hQ0E0TE0T = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', K * numQuad1D, 1);
%   hJmpQ0E0T = hQ0E0Tint - hQ0E0TE0T;
%   normJumpH = sum(reshape((hJmpQ0E0T.^2), numQuad1D, K).' * W(:));
% end

% filename = [problemData.outputBasename '_norms.dat'];
% norm_file = fopen(filename, 'at');
% fprintf(norm_file, '  %16.10e  %16.10e  %16.10e  %16.10e  %16.10e\n', problemData.tEnd, normH, normQ1, normQ2, normJumpH);
% fclose(norm_file);

%% Save final state.
t = problemData.tEnd;
hotstartFile = [ problemData.outputBasename '.mat' ];
save(hotstartFile, 't', 'hDisc', 'q1Disc', 'q2Disc');
fprintf('Saved hotstart data at t=%g to "%s"\n', t, hotstartFile);
end % function

