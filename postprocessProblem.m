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
problemData.t = problemData.tEnd;
problemData = execin([problemData.problemName filesep 'preprocessSubStep'], problemData, problemData.numSteps + 1, 1);

problemData.cDiscRK{end, 3} = reshape( problemData.globHQup \ ( ...
                                   problemData.globJu{1} + problemData.globJu{2} + problemData.globVeeJu + problemData.globVeeJh + ...
                                   0.5 * (problemData.globVeeJuRiem + problemData.globVeeJhRiem + problemData.globJuhRiem) + ...
                                   problemData.globJuCoupling{1} + problemData.globJwCoupling + problemData.globKh + ...
                                  (problemData.globHQavg + problemData.globVeeP + problemData.globVeePbdr) * reshape(problemData.cDiscRK{end, 2}.', [], 1) ), ...
                              problemData.N, problemData.g.numT ).';

%% Evaluate errors (if analytical solution given).
if all(isfield(problemData, { 'hCont', 'u1Cont', 'u2Cont' }))
  t = problemData.tEnd;
  hCont = problemData.hCont;
  u1Cont = problemData.u1Cont;
  u2Cont = problemData.u2Cont;
  
  problemData.error = [ computeL2Error1D(problemData.g.g1D, problemData.cDiscRK{end, 1}, ...
                            @(x1) hCont(t, x1), problemData.qOrd + 1, problemData.basesOnQuad1D), ...
                        computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 2}, ...
                            @(x1,x2) u1Cont(t, x1, x2), problemData.qOrd + 1, problemData.basesOnQuad2D), ...
                        computeL2ErrorTetra(problemData.g, problemData.cDiscRK{end, 3}, ...
                            @(x1,x2) u2Cont(t, x1, x2), problemData.qOrd + 1, problemData.basesOnQuad2D) ];

  fprintf('L2 errors of h, u1, u2 w.r.t. the analytical solution: %g, %g, %g\n', problemData.error);
end % if

%% Visualize final state.
if problemData.isVisGrid, visualizeGridTetra(problemData.g); end

if problemData.isVisSol
  nOutput = ceil(problemData.numSteps / problemData.outputFrequency);
  cLagr = { projectDataDisc2DataLagrTensorProduct(problemData.cDiscRK{end, 2}), ...
            projectDataDisc2DataLagrTensorProduct(problemData.cDiscRK{end, 3}) };
  visualizeDataLagrTetra(problemData.g, cLagr, {'u1', 'u2'}, problemData.outputBasename, nOutput, problemData.outputTypes, struct('velocity', {{'u1','u2'}}));
end % if

%% Save final state.
t = problemData.tEnd;
hDisc = problemData.cDiscRK{end, 1}; %#ok<NASGU>
u1Disc = problemData.cDiscRK{end, 2}; %#ok<NASGU>
u2Disc = problemData.cDiscRK{end, 3}; %#ok<NASGU>
hotstartFile = [ problemData.outputBasename '.mat' ];
save(hotstartFile, 't', 'hDisc', 'u1Disc', 'u2Disc');
fprintf('Saved hotstart data at t=%g to "%s"\n', t, hotstartFile);
end % function

