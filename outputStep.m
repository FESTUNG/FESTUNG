% Last step of the four-part algorithm in the main loop.

%===============================================================================
%> @file sweVert/outputStep.m
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
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed last in each loop iteration and writes output
%> files that can later be visualized using TecPlot, Paraview, or others,
%> depending on the chosen file types in configureProblem().
%>
%> @param  problemData  A struct with problem parameters, precomputed
%>                      fields, and solution data structures, as provided 
%>                      by configureProblem() and preprocessProblem(). 
%>                      @f$[\text{struct}]@f$
%> @param  nStep        The current iteration number of the main loop. 
%>
%> @retval problemData  The input struct enriched with post-processed data
%>                      for this loop iteration. @f$[\text{struct}]@f$
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
    visualizeDataLagrTetra(problemData.g, cLagr(2:3), {'u1', 'u2'}, problemData.outputBasename, nOutput, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
    isOutput = true;
  end % if
  
  if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
    t = problemData.t(1);
    hCont = problemData.hCont;
    u1Cont = problemData.u1Cont;
    u2Cont = problemData.u2Cont;
        
    if problemData.isVisSol
      problemData.pdExact.cDiscRK{end, 1} = projectFuncCont2DataDisc1D(problemData.pdExact.g.g1D, @(x1) hCont(t, x1), problemData.qOrd, problemData.barHatM, problemData.basesOnQuad1D);

      problemData.pdExact = problemData.fn_adaptFreeSurface(problemData.pdExact);

      problemData.pdExact.cDiscRK{end, 2} = projectFuncCont2DataDiscTetra(problemData.pdExact.g, @(x1,x2) u1Cont(t, x1, x2), problemData.qOrd, problemData.pdExact.globM, problemData.basesOnQuad2D);
      problemData.pdExact.cDiscRK{end, 3} = projectFuncCont2DataDiscTetra(problemData.pdExact.g, @(x1,x2) u2Cont(t, x1, x2), problemData.qOrd, problemData.pdExact.globM, problemData.basesOnQuad2D);

      cLagr = { projectDataDisc2DataLagrTensorProduct(problemData.pdExact.cDiscRK{end, 2}), ...
                projectDataDisc2DataLagrTensorProduct(problemData.pdExact.cDiscRK{end, 3}) };

      visualizeDataLagrTetra(problemData.pdExact.g, cLagr, {'u1', 'u2'}, [ problemData.outputBasename '_ex' ], nOutput, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
    end % if
    
    problemData.error = [ computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{1, 1}, ...
                              @(x1) hCont(t, x1), problemData.qOrd + 1, problemData.basesOnQuad1D), ...
                          computeL2ErrorTetra(problemData.g, problemData.cDiscRK{1, 2}, ...
                              @(x1,x2) u1Cont(t, x1, x2), problemData.qOrd + 1, problemData.basesOnQuad2D), ...
                          computeL2ErrorTetra(problemData.g, problemData.cDiscRK{2, 3}, ...
                              @(x1,x2) u2Cont(t, x1, x2), problemData.qOrd + 1, problemData.basesOnQuad2D) ];

    fprintf('L2 errors of h, u1, u2 w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
    isOutput = true;
  end % if
  
  if ~isOutput && nStep > problemData.outputFrequency
    fprintf(repmat('\b', 1, 11));
  end % if
  fprintf('%3.0f %% done\n', nStep / problemData.numSteps * 100);
end % if
end % function