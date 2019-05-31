% Last step of the four-part algorithm in the main loop.

%===============================================================================
%> @file
%>
%> @brief Last step of the four-part algorithm in the main loop.
%===============================================================================
%>
%> @brief Last step of the four-part algorithm in the main loop.
%>
%> The main loop repeatedly executes four steps until the number of
%> iterations provided by configureProblem in the parameter
%> <code>numSteps</code> is reached. These four steps are:
%>
%>  1. preprocessStep()
%>  2. solveStep()
%>  3. postprocessStep()
%>  4. outputStep()
%> 
%> This routine is executed last in each loop iteration and is intended to
%> provide output operations for the solution, e.g., write it to a file
%> for later visualization.
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
%> @copyright 2014-2016 Hennes Hajduk, Balthasar Reuter, Florian Frank, Vadym Aizinger
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
%% visualization
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species} && (mod(nStep, problemData.outputFrequency{species}) == 0 || problemData.isFinished)
    cLagrange = projectDataDisc2DataLagr(problemData.concDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                      problemData.outputBasename{species}, ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
    if problemData.isMask(species)
      visualizeDataLagr(problemData.g, problemData.mask(:,species), ['mask_' num2str(species)], ...
                        [problemData.outputBasename{species} '_mask'], ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
    end % if
  end % if
end % for
% TODO use this kind of visualization
% isVis = cell2mat( cellfun(@(vis, freq) vis && (mod(nStep, freq) == 0 || problemData.isFinished), ...
%                            problemData.isVisSol, problemData.outputFrequency, 'UniformOutput', false) );
% if any(isVis)
%   isVisMask = problemData.isMask(isVis);
%   
%   [~, varName] = cellfun(@(str) fileparts(str), problemData.outputBasename(isVis), 'UniformOutput', false);
%   varName = [ varName, cellfun(@(str) [str '_mask'], varName(isVisMask), 'UniformOutput', false) ]';
%   
%   cLagr = [ cellfun(@(c) projectDataDisc2DataLagr(c), problemData.concDisc(isVis), 'UniformOutput', false); ...
%             mat2cell(problemData.mask(:, isVis & problemData.isMask), size(problemData.mask, 1), ones(size(isVis & problemData.isMask)))' ];
%   
%   visualizeDataLagr(problemData.g, cLagr, varName, problemData.outputBasename{1}, nStep, problemData.outputTypes{1});
% end % if
%% Mass balance
if problemData.isCheckMass
  mass = cellfun(@(c) sum(c(:,1)), problemData.concDisc, 'UniformOutput', false);
  dlmwrite('mass.csv', [nStep, cell2mat(mass)'], '-append', 'Precision', 18)
  dlmwrite('mass_bdr.csv', [nStep, problemData.massLossBdr], '-append', 'Precision', 18)
end % if
end % function
