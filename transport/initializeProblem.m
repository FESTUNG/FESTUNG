% Fills the problem's data structures with initial data (if applicable).

%===============================================================================
%> @file template/initializeProblem.m
%>
%> @brief Fills the problem's data structures with initial data.
%===============================================================================
%>
%> @brief Fills the problem's data structures with initial data.
%>
%> This routine is called after template/preprocessProblem.m.
%>
%> Before entering the main loop of the solution algorithm, this routine
%> fills the problem's data structures with initial data.
%>
%> It loads hotstart data from a file or projects initial data for the
%> primary variables.
%> The initial state of the system is visualized.
%>
%> @param  problemData  A struct with problem parameters and precomputed
%>                      fields, as provided by configureProblem() and 
%>                      preprocessProblem(). @f$[\text{struct}]@f$
%>
%> @retval problemData  The input struct enriched with initial data.
%>                      @f$[\text{struct}]@f$
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
function problemData = initializeProblem(problemData)
%% Initial data.
problemData.cDisc = cell(problemData.numSpecies,1);

if problemData.isMask
  problemData.mask = false(problemData.K, problemData.numSpecies);
else
  problemData.mask = true(problemData.K, problemData.numSpecies);
end % if

problemData.numOperations = 0;

hDisc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.hCont(0,x1,x2), 2*problemData.p+1,problemData.hatM, problemData.basesOnQuad);
hQ0T = (hDisc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.');

for species = 1:problemData.numSpecies
  problemData.cDisc{species} = projectFuncCont2DataDisc(problemData.g, problemData.cH0Cont{species}, 2*problemData.p+1, ...
                                                        problemData.hatM, problemData.basesOnQuad);
  if problemData.isSlopeLim{species}
    cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont{species}(0, x1, x2));
    [problemData.cDisc{species}, minMaxV0T] = applySlopeLimiterDisc(problemData.g, problemData.cDisc{species}, ...
                                                                    problemData.g.markV0TbdrD, cDV0T, problemData.globM, ...
                                                                    problemData.globMDiscTaylor, problemData.basesOnQuad, ...
                                                                    problemData.typeSlopeLim{species});
    if problemData.isMask
      problemData.mask(:,species) = (max(minMaxV0T{2},[],2) - min(minMaxV0T{1},[],2)) >= problemData.maskTol;
      
    end % if
  end % if
  
  % Initial error TODO maybe for non-integrated concentrations too?
  fprintf('L2 error w.r.t. the initial condition: %g\n', ...
    computeL2Error(problemData.g, problemData.cDisc{species}, problemData.cH0Cont{species}, 2*problemData.p, problemData.basesOnQuad));
  
  % Visualization of inital condition.
  if problemData.isVisSol{species}
    cLagrange = projectDataDisc2DataLagr(problemData.cDisc{species});
    visualizeDataLagr(problemData.g, cLagrange, ['cH_' num2str(species) '_h'], ... 
                      [problemData.outputBasename{species} '_H'], 0, problemData.outputTypes{species})
    
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.') ./ hQ0T;
    dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
    cLagrange = projectDataDisc2DataLagr(dataDisc);
    visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                      problemData.outputBasename{species}, 0, problemData.outputTypes{species});
    if problemData.isMask
      visualizeDataLagr(problemData.g, problemData.mask(:,species), ['mask_' num2str(species)], ...
                        [problemData.outputBasename{species} '_mask'], 0, problemData.outputTypes{species});
    end % if
  end % if
end % for
%% Initialize time stepping.
problemData.isFinished = false;
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', ...
  problemData.tEnd, problemData.tau, problemData.numSteps)
end % function