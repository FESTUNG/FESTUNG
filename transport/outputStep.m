% Last step of the four-part algorithm in the main loop.

%===============================================================================
%> @file template/outputStep.m
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
hQ0T = (problemData.hDisc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.');
for species = 1:problemData.numSpecies
  if problemData.isVisSol{species} && mod(nStep, problemData.outputFrequency{species}) == 0
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['cH_' num2str(species) '_h'], ...
                      [problemData.outputBasename{species} '_H'], ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.') ./ hQ0T;
    dataDisc = problemData.swe_projectDataQ0T2DataDisc(dataQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
    cLagrange = projectDataDisc2DataLagr(dataDisc);
    visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                      problemData.outputBasename{species}, ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
    if problemData.isMask
      visualizeDataLagr(problemData.g, problemData.mask(:,species), ['mask_' num2str(species)], ...
                        [problemData.outputBasename{species} '_mask'], ceil(nStep / problemData.outputFrequency{species}), problemData.outputTypes{species});
    end % if
  end % if
end % for
if ~isempty(dataLagr)
  visualizeDataLagr(problemData.g, dataLagr, varName, problemData.outputBasename, nStep, problemData.outputTypes)
end % if
end % function
