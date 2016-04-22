% Fills the problem's data structures with initial data (if applicable).

%===============================================================================
%> @file template/initializeProblem.m
%>
%> @brief Fills the problem's data structures with initial data (if applicable).
%===============================================================================
%>
%> @brief Fills the problem's data structures with initial data (if applicable).
%>
%> This routine is called after template/preprocessProblem.m.
%>
%> Before entering the main loop of the solution algorithm, this routine
%> fills the problem's data structures with initial data.
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
%> @copyright 2014-2016 Balthasar Reuter, Florian Frank, Vadym Aizinger
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
K = problemData.K;
N = problemData.N;

%% Primary unknowns H, uH, vH (each K*N).
if problemData.isSolutionAvail
  h0Cont = @(x1,x2) problemData.xiCont(x1,x2,problemData.t0) - problemData.zbCont(x1,x2);
  uH0Cont = @(x1,x2) h0Cont(x1,x2) .* problemData.uCont(x1,x2,problemData.t0);
  vH0Cont = @(x1,x2) h0Cont(x1,x2) .* problemData.vCont(x1,x2,problemData.t0);
else
  h0Cont = @(x1,x2)  -problemData.zbCont(x1,x2);
  uH0Cont = @(x1,x2) zeros(size(x1));
  vH0Cont = @(x1,x2) zeros(size(x1));
end % if
  
problemData.cDisc = zeros(K,N,3);
problemData.cDisc(:,:,1) = projectFuncCont2DataDisc(problemData.g, h0Cont, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
problemData.cDisc(:,:,2) = projectFuncCont2DataDisc(problemData.g, uH0Cont, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);
problemData.cDisc(:,:,3) = projectFuncCont2DataDisc(problemData.g, vH0Cont, 2*p, problemData.refElemPhiPhi, problemData.basesOnQuad);

problemData.cDisc(:,:,1) = correctMinValueExceedanceDisc(problemData.cDisc(:,:,1), problemData.sysMinValueCorrection, 0, problemData.minValueHeight, 1000);

%% Visualize initial solution.
visualizeSolution(problemData, 0);

%% Initialize waitbar.
if problemData.isWaitbar
  str  = strcat( ['% done. Simulating refinement level =', ' ', num2str(problemData.refinement), ', p =', ' ', num2str(p), '.' ]);
  problemData.waitbar = waitbar(0, strcat([ 'Time stepping:', ' ', num2str(0), str ]));
end % if
end % function