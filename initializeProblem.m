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

problemData.numOperations = 0;

if ~isfield(problemData, 'h0Disc') % TODO
  problemData.h0Disc = projectFuncCont2DataDisc(problemData.g, @(x1,x2) problemData.hCont(0,x1,x2), 2*problemData.p+1,problemData.hatM, problemData.basesOnQuad);
end % if
hQ0T = problemData.h0Disc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.';

for species = 1:problemData.numSpecies
  if ~isfield(problemData, 'cH0Disc') % TODO
    problemData.cDisc{species} = projectFuncCont2DataDisc(problemData.g, problemData.cH0Cont{species}, 2*problemData.p+1, ...
                                                          problemData.hatM, problemData.basesOnQuad);
  else
    problemData.cDisc{species} = problemData.cH0Disc{species};
  end % if
  
  if problemData.isVisSol{species} || problemData.isSlopeLim{species}
    
    % Compute the concentration
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.') ./ hQ0T;
    dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
  
    % Limiting the concentration
    if problemData.isSlopeLim{species}

      cDV0T = computeFuncContV0T(problemData.g, @(x1, x2) problemData.cDCont{species}(0, x1, x2));
      [dataDisc, minMaxV0T] = applySlopeLimiterDisc(problemData.g, dataDisc, problemData.g.markV0TbdrD, cDV0T, problemData.globM, ...
                                                    problemData.globMDiscTaylor, problemData.basesOnQuad, problemData.typeSlopeLim{species});
      
      % Compute the integrated concentration
      dataDiscQ0T = dataDisc * problemData.basesOnQuad.phi2D{max(2*problemData.p,1)}.';
      problemData.cDisc{species} = projectDataQ0T2DataDisc(dataDiscQ0T .* hQ0T, 2*problemData.p, problemData.hatM, problemData.basesOnQuad);
      
      if problemData.isMask(species)
        problemData.mask(:,species) = computeMask(minMaxV0T, problemData.maskTol, problemData.maskType);
        
        if isequal(problemData.mask(:,species), zeros(problemData.K,1)) % possible workaround
          problemData.mask(1,species) = true;
        end % if
      else
        problemData.mask(:,species) = true(problemData.K, 1);
      end % if
    else
      problemData.mask(:,species) = true(problemData.K, 1);
    end % if
    
    if problemData.isVisSol{species}
      cLagrange = projectDataDisc2DataLagr(dataDisc);
      visualizeDataLagr(problemData.g, cLagrange, ['c_' num2str(species) '_h'], ...
                        problemData.outputBasename{species}, 0, problemData.outputTypes{species});

      if problemData.isMask(species)
        if isequal(problemData.mask(:,species), [1; zeros(problemData.K-1,1)]) % TODO due to the workaround
          problemData.mask(1,species) = 0;
        end % if
        visualizeDataLagr(problemData.g, problemData.mask(:,species), ['mask_' num2str(species)], ...
                          [problemData.outputBasename{species} '_mask'], 0, problemData.outputTypes{species});
        if isequal(problemData.mask(:,species), zeros(problemData.K,1)) % TODO due to the workaround
          problemData.mask(1,species) = 1;
        end % if
      end % if
    end % if
  else
    problemData.mask(:,species) = true(problemData.K, 1);
  end % if
end % for
%% Initialize time stepping.
problemData.isFinished = false;
fprintf('Starting time integration from 0 to %g using time step size %g (%d steps).\n', problemData.tEnd, problemData.tau, problemData.numSteps)
end % function