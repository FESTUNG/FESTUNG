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
%>  1. @link swe_2dv/preprocessStep.m @endlink
%>  2. @link swe_2dv/solveStep.m @endlink
%>  3. @link swe_2dv/postprocessStep.m @endlink
%>  4. @link swe_2dv/outputStep.m @endlink
%> 
%> This routine is executed last in each loop iteration and writes output
%> files that can later be visualized using TecPlot, Paraview, or others,
%> depending on the chosen file types in @link swe_2dv/configureProblem.m @endlink.
%>
%> If analytical solution data for @f$h, u^1, u^2@f$ are given, the L2-errors
%> of the current state are evaluated and printed.
%>
%> Furthermore, the current process with respect to the number of executed time
%> steps is printed.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields (either filled with initial data or the solution
%>                      from the previous loop iteration), as provided by 
%>                      @link swe_2dv/configureProblem.m @endlink and 
%>                      @link swe_2dv/preprocessProblem.m @endlink. 
%%>                     @f$[\text{struct}]@f$
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
%% Visualization
if mod(nStep-1, problemData.outputFrequency) == 0
  isOutput = false;
  nOutput = (nStep-1) / problemData.outputFrequency;
  
  if problemData.isVisSol
    cLagr = { projectDataDisc2DataLagr1D(problemData.cDiscRK{1, 1}), ...
              projectDataDisc2DataLagrTensorProduct(problemData.cDiscRK{1, 2}), ...
              projectDataDisc2DataLagrTensorProduct(problemData.cDiscRK{2, 3}) };
    visualizeDataLagr1D(problemData.g.g1D, cLagr{1}, { 'h' }, [ problemData.outputBasename '_h' ], nOutput, problemData.outputTypes);
    visualizeDataLagrQuadri(problemData.g, cLagr(2:3), {'u1', 'u2'}, problemData.outputBasename, nOutput, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
    isOutput = true;
  end % if
  
  if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
    t = problemData.t(1);
    hCont = problemData.hCont;
    u1Cont = problemData.u1Cont;
    u2Cont = problemData.u2Cont;
        
    if problemData.isVisSol
      problemData.pdExact.cDiscRK{end, 1} = projectFuncCont2DataDisc1D(problemData.pdExact.g.g1D, @(x1) hCont(t, x1), problemData.qOrd, problemData.hatBarM, problemData.basesOnQuad1D);

      problemData.pdExact = problemData.fn_adaptFreeSurface(problemData.pdExact);

      problemData.pdExact.cDiscRK{end, 2} = projectFuncCont2DataDiscQuadri(problemData.pdExact.g, @(x1,x2) u1Cont(t, x1, x2), problemData.qOrd, problemData.pdExact.globM, problemData.basesOnQuad2D);
      problemData.pdExact.cDiscRK{end, 3} = projectFuncCont2DataDiscQuadri(problemData.pdExact.g, @(x1,x2) u2Cont(t, x1, x2), problemData.qOrd, problemData.pdExact.globM, problemData.basesOnQuad2D);

      cLagr = { projectDataDisc2DataLagrTensorProduct(problemData.pdExact.cDiscRK{end, 2}), ...
                projectDataDisc2DataLagrTensorProduct(problemData.pdExact.cDiscRK{end, 3}) };

      visualizeDataLagrQuadri(problemData.pdExact.g, cLagr, {'u1', 'u2'}, [ problemData.outputBasename '_ex' ], nOutput, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
    end % if
    
    problemData.error = [ computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{1, 1}, ...
                              @(x1) hCont(t, x1), problemData.qOrd + 1, problemData.basesOnQuad1D), ...
                          computeL2ErrorQuadri(problemData.g, problemData.cDiscRK{1, 2}, ...
                              @(x1,x2) u1Cont(t, x1, x2), problemData.qOrd + 1, problemData.basesOnQuad2D), ...
                          computeL2ErrorQuadri(problemData.g, problemData.cDiscRK{2, 3}, ...
                              @(x1,x2) u2Cont(t, x1, x2), problemData.qOrd + 1, problemData.basesOnQuad2D) ];

    fprintf('L2 errors of h, u1, u2 w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
    isOutput = true;
  end % if
  
  if ~isOutput && nStep > problemData.outputFrequency
    fprintf(repmat('\b', 1, 11));
  end % if
  fprintf('%3.0f %% done\n', nStep / problemData.numSteps * 100);
end % if

%% Norms for stability estimate
% if mod(nStep-1, 100) == 0
%   K = problemData.g.numT;
%   [Q, W] = quadRule1D(problemData.qOrd); numQuad1D = length(Q);
  
%   xiDisc = reshape((problemData.cDiscRK{1, 1} + problemData.zBotDisc).', [], 1);
%   normXi = xiDisc.' * (problemData.globBarM * xiDisc);
  
%   uDisc = reshape(problemData.cDiscRK{1, 2}.', [], 1);
%   normU = uDisc.' * (problemData.globM * uDisc);
  
%   u1Q0E0Tint = reshape(problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, 3) * bsxfun( ...
%                       @times, problemData.g.markE0Tint(:, 3), problemData.cDiscRK{1, 2}).', ...
%                   K * numQuad1D, 1);
%   cDiscThetaPhi = problemData.basesOnQuad2D.phi1D{problemData.qOrd}(:, :, mapLocalEdgeIndexQuadri(3)) * problemData.cDiscRK{1, 2}.';
%   u1Q0E0TE0T = reshape(cDiscThetaPhi * problemData.g.markE0TE0T{3}.', K * numQuad1D, 1);
%   u1JmpQ0E0T = u1Q0E0Tint - u1Q0E0TE0T;
%   normJumpU = sum(problemData.g.areaE0T{3} .* (reshape((u1JmpQ0E0T.^2), numQuad1D, K).' * W(:)));
  
%   filename = [problemData.outputBasename '_norms.dat'];
%   norm_file = fopen(filename, 'at');
%   fprintf(norm_file, '  %16.10e  %16.10e  %16.10e  %16.10e  %16.10e  %16.10e\n', problemData.t(1), normXi, normU, problemData.normQ(1), problemData.normQ(2), normJumpU);
%   fclose(norm_file);
% end % if
end % function