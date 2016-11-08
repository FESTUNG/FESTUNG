% Performs all post-processing steps, such as error estimates of the final
% solution, etc.

%===============================================================================
%> @file template/postprocessProblem.m
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%===============================================================================
%>
%> @brief Performs all post-processing tasks, such as error estimates of the 
%>        final solution, etc.
%>
%> This routine is called after the main loop.
%>
%> It can include things such as error estimates, an output operation of
%> the final solution, etc.
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
function problemData = postprocessProblem(problemData)
p = problemData.p;
qOrd2D = max(2*p,1);

if problemData.isSolutionAvailable
  problemData.errors = zeros(problemData.numSpecies,2);
  for species = 1 : problemData.numSpecies
    %% Error evaluation
    problemData.errors(species,1) = computeL2Error(problemData.g, problemData.cDisc{species}, @(x1,x2) problemData.solCont{species}(problemData.tEnd,x1,x2) .* problemData.hCont(problemData.tEnd,x1,x2), 2*p, problemData.basesOnQuad);
    dataQ0T = (problemData.cDisc{species} * problemData.basesOnQuad.phi2D{qOrd2D}.') ./ problemData.hQ0T;
    dataDisc = projectDataQ0T2DataDisc(dataQ0T, 2*p, problemData.hatM, problemData.basesOnQuad);
    problemData.errors(species,2) = computeL2Error(problemData.g, dataDisc, @(x1,x2) problemData.solCont{species}(problemData.tEnd,x1,x2), 2*p, problemData.basesOnQuad);
  end % for
end % if

if ~isequal(problemData.isMask, zeros(problemData.numSpecies, 1))
  fprintf([ '%d operations were needed for solving the transport model on %d triangles, %d time steps and %d species.\n' ...
            'This corresponds to %3.1f %% of the full number of operations.\n'], ...
            problemData.numOperations, problemData.g.numT, problemData.numSteps, problemData.numSpecies, ...
            problemData.numOperations / (problemData.g.numT*problemData.numSteps*problemData.numSpecies) * 100);
end % if
end % function

