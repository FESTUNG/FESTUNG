% Last step of the four-part algorithm in the main loop.

%===============================================================================
%> @file
%>
%> @brief Last step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Last step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the parameter
%> <code>problemData.isFinished</code> becomes <code>true</code>.
%> These four steps are:
%>
%>  1. darcy_2dv/preprocessStep.m
%>  2. darcy_2dv/solveStep.m
%>  3. darcy_2dv/postprocessStep.m
%>  4. darcy_2dv/outputStep.m
%> 
%> This routine is executed last in each loop iteration and writes output
%> files that can later be visualized using TecPlot, Paraview, or others,
%> depending on the chosen file types in darcy_2dv/configureProblem.m.
%>
%> If analytical solution data for @f$h, \vec{q}@f$ are given, the L2-error
%> of the current state is evaluated and printed.
%>
%> Furthermore, the current process with respect to the number of executed time
%> steps is printed.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      darcy_2dv/configureProblem.m and 
%>                      darcy_2dv/preprocessProblem.m. @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with post-processed data
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
function problemData = outputStep(problemData, nStep)
K = problemData.g.numT;
N = problemData.N;
%% Visualization
if mod(nStep, problemData.outputFrequency) == 0
  isOutput = false;
  nOutput = nStep / problemData.outputFrequency;
  
  if problemData.isVisSol
    cLagr = { projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)'), ...
              projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(1 : K*N), N, K)'), ...
              projectDataDisc2DataLagrTensorProduct(reshape(problemData.sysY(K*N+1 : 2*K*N), N, K)') };
    visualizeDataLagrQuadri(problemData.g, cLagr, {'h', 'q1', 'q2'}, problemData.outputBasename, nOutput, problemData.outputTypes, struct('q', {{'q1','q2'}}));
    isOutput = true;
  end % if
    
  if all(isfield(problemData, { 'hCont', 'q1Cont', 'q2Cont' }))
    t = nStep * problemData.tau;  
    hCont = @(x1,x2) problemData.hCont(t, x1, x2);
    q1Cont = @(x1,x2) problemData.q1Cont(t, x1, x2);
    q2Cont = @(x1,x2) problemData.q2Cont(t, x1, x2);

    problemData.error = [ computeL2ErrorQuadri(problemData.g, reshape(problemData.sysY(2*K*N+1 : 3*K*N), N, K)', ...
                              hCont, problemData.qOrd + 1, problemData.basesOnQuad), ...
                          computeL2ErrorQuadri(problemData.g, reshape(problemData.sysY(1 : K*N), N, K)', ...
                              q1Cont, problemData.qOrd + 1, problemData.basesOnQuad), ...
                          computeL2ErrorQuadri(problemData.g, reshape(problemData.sysY(K*N+1 : 2*K*N), N, K)', ...
                              q2Cont, problemData.qOrd + 1, problemData.basesOnQuad) ];

    fprintf('L2 errors of h, q1, q2 w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
    
    if problemData.isVisSol
      cDisc = { projectFuncCont2DataDiscQuadri(problemData.g, hCont, problemData.qOrd, problemData.globM, problemData.basesOnQuad), ...
                projectFuncCont2DataDiscQuadri(problemData.g, q1Cont, problemData.qOrd, problemData.globM, problemData.basesOnQuad), ...
                projectFuncCont2DataDiscQuadri(problemData.g, q2Cont, problemData.qOrd, problemData.globM, problemData.basesOnQuad) };
      cLagr = cellfun(@(c) projectDataDisc2DataLagrTensorProduct(c), cDisc, 'UniformOutput', false);
      visualizeDataLagrQuadri(problemData.g, cLagr, {'h', 'q1', 'q2'}, [ problemData.outputBasename '_ex' ], nOutput, problemData.outputTypes, struct('q', {{'q1','q2'}}));
    end % if
    
    isOutput = true;
  end % if
  
  if ~isOutput && nStep > problemData.outputFrequency
    fprintf(repmat('\b', 1, 11));
  end % if
  
  fprintf('%3.0f %% done\n', nStep / problemData.numSteps * 100);
end % if

%% Norms for stability estimate
% if mod(nStep, 10) == 0
%   [Q, W] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);

%   hDisc = problemData.sysY(2*K*N+1 : 3*K*N);
%   normH = hDisc.' * (problemData.globM * hDisc);

%   q1Disc = problemData.sysY(1 : K*N);
%   normQ1 = q1Disc.' * (problemData.globM * q1Disc);

%   q2Disc = problemData.sysY(K*N+1 : 2*K*N);
%   normQ2 = q2Disc.' * (problemData.globM * q2Disc);

%   normJumpH = 0;
%   hDisc = reshape(hDisc, N, K)';
%   for n = 2 : 3
%     hQ0E0Tint = reshape(problemData.basesOnQuad.phi1D{problemData.qOrd}(:, :, n) * bsxfun( ...
%                         @times, problemData.g.markE0Tint(:, n), hDisc).', ...
%                     K * numQuad1D, 1);
%     cDiscThetaPhi = problemData.basesOnQuad.phi1D{problemData.qOrd}(:, :, mapLocalEdgeIndexQuadri(3)) * hDisc.';
%     hQ0E0TE0T = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{n}.', K * numQuad1D, 1);
%     hJmpQ0E0T = hQ0E0Tint - hQ0E0TE0T;
%     normJumpH = sum(reshape((hJmpQ0E0T.^2), numQuad1D, K).' * W(:));
%   end

%   filename = [problemData.outputBasename '_norms.dat'];
%   norm_file = fopen(filename, 'at');
%   fprintf(norm_file, '  %16.10e  %16.10e  %16.10e  %16.10e  %16.10e\n', nStep * problemData.tau, normH, normQ1, normQ2, normJumpH);
%   fclose(norm_file);
% end % if
end % function
