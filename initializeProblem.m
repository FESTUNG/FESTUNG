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
p = problemData.p;
qOrd2D = max(2*p,1);

%% Initial data.
problemData.cDisc = cell(problemData.numSpecies,1);
problemData.cQ0T = cell(problemData.numSpecies,1);
problemData.concDisc = cell(problemData.numSpecies,1);

problemData.numOperations = zeros(1, problemData.numSpecies);
problemData.numElem = zeros(1, problemData.numSpecies);
problemData.minMaxV0T = cell(1, problemData.numSpecies);

if ~isfield(problemData, 'h0Disc')
  problemData.h0Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.hCont(0,x1,x2), 2*p+1,problemData.hatM, problemData.basesOnQuad);
end % if
hQ0T = problemData.h0Disc * problemData.basesOnQuad.phi2D{qOrd2D}.';

for species = 1:problemData.numSpecies
  if ~isfield(problemData, 'cH0Disc')
    problemData.cDisc{species} = projectFuncCont2DataDisc(problemData.g, problemData.cH0Cont{species}, 2*p+1, ...
                                                          problemData.hatM, problemData.basesOnQuad);
  else
    problemData.cDisc{species} = problemData.cH0Disc{species};
  end % if
  
  if problemData.isVisSol{species} || problemData.isSlopeLim{species}
    
    % Compute the concentration
    problemData.cQ0T{species} = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{qOrd2D}.') ./ hQ0T;
    problemData.concDisc{species} = projectDataQ0T2DataDisc(problemData.cQ0T{species}, 2*p, problemData.hatM, problemData.basesOnQuad);
  
    % Limiting the concentration
    if problemData.isSlopeLim{species}

      cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont{species}(0, x1, x2));
      [problemData.concDisc{species}, minMaxV0T] = applySlopeLimiterDisc(problemData.g, problemData.concDisc{species}, ...
                                                                                    problemData.g.markV0TbdrD, cDV0T, problemData.globM, ...
                                                                                    problemData.globMDiscTaylor, problemData.basesOnQuad, ...
                                                                                    problemData.typeSlopeLim{species});
      
      % Compute the integrated concentration
      dataDiscQ0T = problemData.concDisc{species} * problemData.basesOnQuad.phi2D{qOrd2D}.';
      problemData.cDisc{species} = projectDataQ0T2DataDisc(dataDiscQ0T .* hQ0T, 2*p, problemData.hatM, problemData.basesOnQuad);
      
      if problemData.isMask(species)
        problemData.mask(:,species) = computeMask(minMaxV0T, problemData.maskTol(species), problemData.maskType);
      else
        problemData.mask(:,species) = true(problemData.K, 1);
      end % if
    else
      problemData.mask(:,species) = true(problemData.K, 1);
    end % if
    
    if problemData.isVisSol{species} % TODO use different kind of visualization
      cLagrange = projectDataDisc2DataLagr(problemData.concDisc{species});
      visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                        problemData.outputBasename{species}, 0, problemData.outputTypes{species});

      if problemData.isMask(species)
        visualizeDataLagr(problemData.g, problemData.mask(:,species), ['mask_' num2str(species)], ...
                          [problemData.outputBasename{species} '_mask'], 0, problemData.outputTypes{species});
      end % if
    end % if
  else
    problemData.mask(:,species) = true(problemData.K, 1);
  end % if
end % for
problemData.numElem = sum(problemData.mask, 1);

%% Mass balance
if problemData.isCheckMass
%   header = {'nStep', arrayfun(@(n) ['species_' num2str(n)], 1 : problemData.numSpecies, 'UniformOutput', false) };
%   dlmwrite('mass.csv', header)
  dlmwrite('mass.csv', [0, cell2mat(cellfun(@(c) sum(c(:,1)), problemData.concDisc, 'UniformOutput', false))'], 'Precision', 18)
  problemData.massLossBdr = zeros(1, problemData.numSpecies);
%   dlmwrite('mass_bdr.csv', header)
  dlmwrite('mass_bdr.csv', [0, problemData.massLossBdr], 'Precision', 18)
end % if

%% Initialize time stepping.
problemData.isFinished = false;
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', problemData.tEnd, problemData.tau, problemData.numSteps)
end % function